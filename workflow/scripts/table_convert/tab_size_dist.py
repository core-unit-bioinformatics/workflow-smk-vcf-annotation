#!/usr/bin/env python3

import argparse as argp
import collections as col
import re

import pathlib as pl

import pandas as pd
import pysam


# See VCF spec. v4.2 / Section 5.4
# REF ALT Meaning
# s t[p[ piece extending to the right of p is joined after t
# s t]p] reverse comp piece extending left of p is joined after t
# s ]p]t piece extending to the left of p is joined before t
# s [p[t reverse comp piece extending right of p is joined before t
CONST_BND_ORIENTATIONS = {
    ("[", "end"): ("p.SEQ>t.SEQ", "t[p["),
    ("]", "end"): ("SEQ.p>t.RC[SEQ]", "t]p]"),
    ("]", "start"): ("SEQ.p>SEQ.t", "]p]t"),
    ("[", "start"): ("p.SEQ>RC[SEQ].t", "[p[t"),
}


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--vcf",
        "--input",
        "-i",
        dest="vcf",
        type=lambda x: pl.Path(x).resolve(strict=True),
        help="Full path to VCF input file."
    )

    parser.add_argument(
        "--out-table",
        "-tab",
        "-t",
        dest="out_table",
        type=lambda x: pl.Path(x).resolve(),
        default=pl.Path(".").joinpath("vcf_table.tsv.gz"),
        help="Full path to tabular VCF output."
    )

    parser.add_argument(
        "--out-size-dist",
        "-sdist",
        "-s",
        dest="out_sdist",
        type=lambda x: pl.Path(x).resolve(),
        default=pl.Path(".").joinpath("vcf_sdist.tsv"),
        help="Full path to size distribution output."
    )

    parser.add_argument(
        "--sample",
        "-smp",
        type=str,
        dest="sample",
        default=None,
        help="Specify sample name to add to output table."
    )

    parser.add_argument(
        "--callset",
        "-call",
        type=str,
        dest="callset",
        default="callset",
        help="Specify callset name to add to output table."
    )


    args = parser.parse_args()

    return args


class Variant:

    __slots__ = (
        "chrom",
        "start",
        "end",
        "vartype",
        "size",
        "name",
        "sample",
        "callset",
        "vargroup",
        "genotype_1",
        "genotype_2",
        "read_depth",
        "call_quality",
        "genotype_quality",
        "ref_allele",
        "ref_allele_repr",
        "ref_allele_depth",
        "ref_allele_freq",
        "alt_allele",
        "alt_allele_repr",
        "alt_allele_depth",
        "alt_allele_freq",
        "chrom2",
        "chrom2_pos",
        "bnd_join_op",  # the operation to perform for breakend joining
        "bnd_join_vcf",  # the definition as stated in the VCF spec
    )

    def __init__(self, sample=None, callset=None, vcf_record=None, var_desc=None):

        # make sure all slots are set to something
        # that fails for any further operation
        [setattr(self, slot, None) for slot in self.__slots__]

        if all(x is None for x in [sample, callset, vcf_record]):
            assert var_desc is not None
            [setattr(self, key, value) for key, value in var_desc.items()]
        else:
            # set general variant descriptors
            self.chrom = vcf_record.chrom
            self.start = vcf_record.start
            self.end = vcf_record.stop

            self.name = vcf_record.id
            self.sample = sample
            self.callset = callset

            for desc, value in var_desc.items():
                setattr(self, desc, value)

            ref_allele, ref_allele_repr, ref_allele_length = self.normalize_ref_allele(vcf_record.ref)
            self.ref_allele = ref_allele
            self.ref_allele_repr = ref_allele_repr
            alt_allele, alt_allele_repr, alt_allele_length, bnd_spec = self.normalize_alt_allele(vcf_record.alts)
            self.alt_allele = alt_allele
            self.alt_allele_repr = alt_allele_repr

            # set BND infos
            if bnd_spec is not None:
                self.vartype = "BND"
                self.vargroup = "SV"
                if self.chrom2 is None:
                    self.chrom2 = bnd_spec["chrom2"]
                else:
                    assert self.chrom2 == bnd_spec["chrom2"], f"{self.name}: {self.chrom2} / {bnd_spec}"
                self.chrom2_pos = bnd_spec["chrom2_pos"]
                self.bnd_join_op = bnd_spec["bnd_join_op"]
                self.bnd_join_vcf = bnd_spec["bnd_join_vcf"]

            # start auto-completion
            self.auto_infer_size(ref_allele_length, alt_allele_length)
            self.auto_infer_vartype(ref_allele_length, alt_allele_length)

        # this is a general check no matter how
        # the variant record was initialized
        self.auto_complete()

        return

    def __repr__(self):

        var_repr = (
            f"\n=== Sample {self.sample}\n"
            f"=== Callset {self.callset}\n"
            f"=== Variant {self.name}\n"
            "-----\n"
        )
        skip = ["name", "callset", "sample"]
        for slot in self.__slots__:
            if slot in skip:
                continue
            var_repr += f"{slot}: {getattr(self, slot)}\n"
        var_repr += "-----\n"
        return var_repr

    def normalize_ref_allele(self, ref_allele):

        ref_allele_length = None
        ref_allele_repr = "0:0:0:0"
        # default of zero: for incomplete refs such as hg38,
        # some callers produce the ref allele N, which
        # interferes with determining the variant type
        # downstream (ref length is None/not set); this
        # may or may not cause some instances of OBO,
        # i.e. this is a decision how to handle the
        # uncertainty
        ref_allele_length = 0
        assert isinstance(ref_allele, str)
        if ref_allele.upper() in ["N", ".", "*"]:
            ref_allele = "UNK"
        else:
            ref_allele = ref_allele.upper()
            ref_allele_repr = self.build_allele_repr(ref_allele)
            ref_allele_length = len(ref_allele)
            ref_allele = "SEQ"
        if ref_allele_length == 0:
            assert ref_allele == "UNK"
        return ref_allele, ref_allele_repr, ref_allele_length

    def normalize_alt_allele(self, alt_allele):

        alt_allele_length = None
        alt_allele_repr = "0:0:0:0"
        bnd_spec = None
        if isinstance(alt_allele, tuple):
            assert len(alt_allele) == 1
            alt_allele = alt_allele[0]
        if "<" in alt_allele:
            # for symbolic ALTs, the length
            # info has to be inferred from
            # some other source. Make sure
            # that operations are incompatible
            # with the length set here by using
            # None (can't combine None and int)
            alt_allele = "SYM"
        elif alt_allele.upper() in [".", "N", "*"]:
            # same reasoning as for symbolic ALTs
            alt_allele = "UNK"
        elif "[" in alt_allele or "]" in alt_allele:
            # NB: the lookahead here is important to avoid
            # matching allele specs or N if present
            chrom2 = re.search("(chr)?[0-9A-Z]+(?=\:)", alt_allele)
            if chrom2 is None:
                pass
            else:
                chrom2 = chrom2.group(0)
            chrom2_coord = re.search("[0-9]+", alt_allele.split(":")[1])
            if chrom2_coord is None:
                raise ValueError(f"Cannot find BND chrom2 coordinates: {alt_allele}")
            chrom2_coord = int(chrom2_coord.group(0))
            if alt_allele[0] in list("[]"):
                join_type, vcf_rep = CONST_BND_ORIENTATIONS[(alt_allele[0], "start")]
            elif alt_allele[-1] in list("[]"):
                join_type, vcf_rep = CONST_BND_ORIENTATIONS[(alt_allele[-1], "end")]
            else:
                join_type = "UNK"
                vcf_rep = "UNK"
            if join_type == "UNK":
                alt_allele = "UNK"
            else:
                alt_allele = "SEQ"
            alt_allele_length = 0
            bnd_spec = {
                "chrom2": chrom2,
                "chrom2_pos": chrom2_coord,
                "bnd_join_op": join_type,
                "bnd_join_vcf": vcf_rep
            }
        else:
            alt_allele = alt_allele.upper()
            alt_allele_length = len(alt_allele)
            alt_allele_repr = self.build_allele_repr(alt_allele)
            alt_allele = "SEQ"
        return alt_allele, alt_allele_repr, alt_allele_length, bnd_spec

    def build_allele_repr(self, allele):
        nuc_counts = col.Counter(allele.upper())
        allele_repr = (
            f"{nuc_counts['A']}:"
            f"{nuc_counts['C']}:"
            f"{nuc_counts['G']}:"
            f"{nuc_counts['T']}"
        )
        return allele_repr

    def auto_infer_size(self, ref_length, alt_length):

        if self.size is not None:
            pass
        elif isinstance(ref_length, int) and isinstance(alt_length, int):
            if ref_length == alt_length:
                self.size = ref_length
            else:
                self.size = abs(ref_length - alt_length)
        else:
            # make sure that no type casting went wrong
            assert ref_length is None or alt_length is None
            self.size = self.end - self.start
            # can the above cause OBO ?
        if self.size == 0 and self.vartype == "BND":
            # a region in the output BED-like TSV table
            # always has a minimal size of at least 1,
            # so adopt that here for breakends
            self.size = 1
        assert self.size is not None and self.size > 0, f"size affert file: {self}"
        return

    def auto_infer_vartype(self, ref_length, alt_length):

        if self.vartype is not None:
            if self.vargroup is None:
                if self.size < 50:
                    assert self.size > 1
                    self.vargroup = "INDEL"
                else:
                    self.vargroup = "SV"
        elif ref_length is None or alt_length is None:
            # should have been set before
            raise ValueError(f"Cannot determine variant type: {self}")
        else:
            if self.size == 1 and ref_length == alt_length:
                self.vartype = "SNV"
                self.vargroup = "SNV"
            elif 1 < self.size < 50 and ref_length == alt_length:
                # trick, could also be a small INV, but these
                # are usually set via SVTYPE by current callers
                self.vartype = "MNV"
                self.vargroup = "SNV"
            elif 1 <= self.size < 50 and ref_length != alt_length:
                self.vargroup = "INDEL"
                if ref_length > alt_length:
                    self.vartype = "DEL"
                else:
                    self.vartype = "INS"
            elif ref_length > alt_length:
                self.vartype = "DEL"
                self.vargroup = "SV"
            elif ref_length < alt_length:
                self.vartype = "INS"
                self.vargroup = "SV"
            else:
                raise ValueError(f"Cannot determine variant type: {self}")
        assert self.vartype is not None and self.vargroup is not None
        return

    def auto_complete(self):

        for slot in self.__slots__:
            if getattr(self, slot) is not None:
                continue
            if slot == "chrom2":
                if self.vartype != "BND":
                    self.chrom2 = "NA"
                    self.chrom2_pos = -1
                    self.bnd_join_op = "NA"
                    self.bnd_join_vcf = "NA"
                    continue
                else:
                    # cuteSV can report a BND w/o any additional info;
                    # quite implausible (see below, similar for sniffles2),
                    # but it is what it is
                    self.chrom2 = "chrN"
                    self.chrom2_pos = -1
                    self.bnd_join_op = "UNK"
                    self.bnd_join_vcf = "UNK"
                    continue
            if slot == "read_depth":
                assert self.ref_allele_depth is not None and self.alt_allele_depth is not None
                self.read_depth = self.ref_allele_depth + self.alt_allele_depth
                continue
            if slot in ["genotype_quality", "call_quality"]:
                setattr(self, slot, -1)
                continue
            if slot == "ref_allele_freq":
                assert self.read_depth is not None and self.ref_allele_depth is not None
                self.ref_allele_freq = round(self.ref_allele_depth / self.read_depth, 5)
                continue
            if slot == "alt_allele_freq":
                assert self.read_depth is not None and self.alt_allele_depth is not None
                self.alt_allele_freq = round(self.alt_allele_depth / self.read_depth, 5)
                continue
            if slot == "chrom2_pos":
                # unclear how sniffles2 can detect a
                # BND w/o knowing where the
                # other break is on the chromosome 2
                # but it is what it is ...
                self.chrom2_pos = -1
                self.bnd_join_op = "UNK"
                self.bnd_join_vcf = "UNK"
                continue
            raise ValueError(f"Unset variant descriptor with no default: {slot} / {self}")
        return

    def to_table_row(self):
        row = tuple(getattr(self, slot) for slot in self.__slots__)
        return row

    @classmethod
    def from_series(cls, init_series):
        if isinstance(init_series, tuple):
            # assume row/series index
            init_series = init_series[1]
        assert isinstance(init_series, pd.Series)
        return Variant.from_dict(init_series.to_dict())

    @classmethod
    def from_dict(cls, init_dict):
        return Variant(sample=None, callset=None, vcf_record=None, var_desc=init_dict)

    @classmethod
    def get_header(cls):
        return cls.__slots__


def parser_genotype(genotype):
    assert len(genotype) == 2
    try:
        gt1 = int(genotype[0])
    except TypeError:
        gt1 = -1
    try:
        gt2 = int(genotype[1])
    except TypeError:
        gt2 = -1
    return gt1, gt2


def parse_int_tuple(value):
    # no clue why the parser
    # returns that as a tuple
    # for pbmm2 callsets

    if isinstance(value, tuple):
        assert len(value) == 1
        value = abs(int(value[0]))
    else:
        value = abs(int(value))
    return value


def parse_float_tuple(value):
    # no clue why the parser
    # returns that as a tuple
    # for pbmm2 callsets

    if isinstance(value, tuple):
        assert len(value) == 1
        value = round(float(value[0]), 5)
    else:
        value = round(float(value), 5)
    return value


def main():

    args = parse_command_line()

    vcf_key_map = {
        "QUAL": ("call_quality", int),
        "GQ": ("genotype_quality", int),
        "DP": ("read_depth", int),
        "DR": ("ref_allele_depth", int),
        "DV": ("alt_allele_depth", int),
        "VAF": ("alt_allele_freq", parse_float_tuple),
        "AF": ("alt_allele_freq", parse_float_tuple),
        "RE": ("alt_allele_depth", int),
        "SUPPORT": ("alt_allele_depth", int),
        "GT": ("genotype", parser_genotype),
        "SVLEN": ("size", parse_int_tuple),
        "SVTYPE": ("vartype", str.upper),
        "CHR2": ("chrom2", str),
        "AD": ("allele_depth", None)  # dummy just to keep the data
    }

    sample = args.sample
    callset = args.callset

    table = []

    with pysam.VariantFile(args.vcf) as vcf:
        for record in vcf:
            vcf_samples = [vcf_sample for vcf_sample in record.samples]
            if len(vcf_samples) != 1:
                raise ValueError("Only single-sample VCF supported")
            vcf_sample = vcf_samples[0]
            # if the user did no specify as sample
            # name, use the one from the VCF
            if sample is None:
                sample = vcf_sample

    with pysam.VariantFile(args.vcf) as vcf:
        format_fields = sorted(vcf.header.formats.keys())
        format_fields = [
            fmt for fmt in format_fields
            if fmt in vcf_key_map
        ]
        info_fields = sorted(vcf.header.info.keys())
        info_fields = [
            info for info in info_fields
            if info in vcf_key_map
        ]
        for record in vcf:
            variant_desc = {}
            try:
                variant_desc["call_quality"] = int(record.qual)
            except TypeError:
                variant_desc["call_quality"] = -1

            for fmt_field, fmt_value in record.samples[vcf_sample].items():
                if fmt_field not in format_fields:
                    continue
                if fmt_field == "AD":
                    assert len(fmt_value) == 2, f"{fmt_value}"
                    ref_allele_depth, var_allele_depth = fmt_value
                    variant_desc["ref_allele_depth"] = int(ref_allele_depth)
                    variant_desc["alt_allele_depth"] = int(var_allele_depth)
                    continue
                if fmt_field == "GT":
                    assert len(fmt_value) == 2
                    gt1, gt2 = parser_genotype(fmt_value)
                    variant_desc["genotype_1"] = gt1
                    variant_desc["genotype_2"] = gt2
                    continue
                key, value_parser = vcf_key_map[fmt_field]
                try:
                    variant_desc[key] = value_parser(fmt_value)
                except TypeError:
                    print(fmt_field, fmt_value)
                    print(record.id)
                    raise
            for info_field, info_value in record.info.items():
                if not (info_field in format_fields or info_field in info_fields):
                    continue
                try:
                    key, value_parser = vcf_key_map[info_field]
                    variant_desc[key] = value_parser(info_value)
                except KeyError:
                    pass
                except TypeError:
                    print(info_field, info_value)
                    print(record.id)
                    raise
            # note here: use the sample name specified by the user
            variant = Variant(sample, callset, record, variant_desc)
            table.append(variant.to_table_row())

    table = pd.DataFrame.from_records(table, columns=Variant.get_header())

    stats = table.groupby(
        ["vargroup", "vartype"]
    )["size"].describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])

    size_agg = table.groupby(
        ["vargroup", "vartype"]
    )["size"].sum()

    stats["mean"] = stats["mean"].round(2)
    stats["std"] = stats["std"].round(2)
    stats = stats.join(size_agg, how="outer")

    rename_columns = {
        "mean": "size_mean",
        "min": "size_min",
        "std": "size_stddev",
        "5%": "size_5pctile",
        "25%": "size_25pctile",
        "50%": "size_median",
        "75%": "size_75pctile",
        "95%": "size_95pctile",
        "max": "size_max",
        "size": "size_total_bp"
    }
    stats.rename(rename_columns, axis=1, inplace=True)

    args.out_table.parent.mkdir(exist_ok=True, parents=True)
    table.to_csv(args.out_table, sep="\t", header=True, index=False)

    args.out_sdist.parent.mkdir(exist_ok=True, parents=True)
    stats.to_csv(args.out_sdist, sep="\t", header=True, index=True)

    return 0


if __name__ == "__main__":
    main()
