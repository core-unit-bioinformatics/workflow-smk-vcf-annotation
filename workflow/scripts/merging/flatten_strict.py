#!/usr/bin/env python3

import argparse as argp
import collections as col
import contextlib as ctl
import csv
import functools as fnt
import hashlib as hl
import io
import gzip
import operator as op
import pathlib as pl

import xopen


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--merged-table",
        "-b",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="merged_table",
        help="Path to merged variants (output of bedtools merge).",
        required=True
    )

    parser.add_argument(
        "--variant-group",
        "-g",
        choices=["INDEL", "SNV"],
        type=str,
        dest="variant_group",
        help="Type of variants in merged input table.",
        required=True
    )

    parser.add_argument(
        "--singletons",
        "-s",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="singletons",
        help="Path to output file for singletons."
    )

    parser.add_argument(
        "--multicalls",
        "-m",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="multicalls",
        help="Path to output file for multicalls / confirmed calls."
    )

    parser.add_argument(
        "--count-stats",
        "-c",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="count_stats",
        help="Path to output file for counting statistics."
    )

    args = parser.parse_args()

    return args


def get_table_header(variant_group):

    positional = [
        "chrom", "start", "end"
    ]

    # columns 6,7,8,16,20
    commons = [
        "name", "sample", "callset_id", "ref_allele", "alt_allele"
    ]

    known_headers = {
        "SNV": positional + commons,
        "INDEL": positional + ["var_type", "var_length"] + commons
    }

    return known_headers[variant_group]


def output_header(output_type):

    positional = [
        "chrom", "start", "end"
    ]

    commons = [
        "name", "var_type", "var_length",
        "sample", "callset_id"
    ]

    only_multicalls = ["group_id", "group_size"]

    allele_rep = ["ref_code", "alt_code"]

    if output_type == "singletons":
        header = positional + commons + allele_rep
    elif output_type == "multicalls":
        header = positional + commons + only_multicalls + allele_rep
    else:
        raise ValueError(f"Unknown output type: {output_type}")
    return header


def check_or_get_column_header(table_file, variant_group):

    first_line = table_file.readline()
    if first_line.startswith("#"):
        check_header = first_line.strip("#").strip().split()
        header = []
        for item in check_header:
            if item == "call_id":
                header.append("name")
            else:
                header.append(item)
    else:
        table_file.seek(0)
        header = get_table_header(variant_group)
    return header


def field_getter(variant_group, column_names):

    if variant_group == "SNV":
        select_fields = ["name", "sample", "callset_id", "ref_allele", "alt_allele"]
    elif variant_group == "INDEL":
        select_fields = ["name", "sample", "callset_id", "var_type", "var_length", "ref_allele", "alt_allele"]
    else:
        raise NotImplementedError(f"no field getter for: {variant_group}")

    if not all(field in column_names for field in column_names):
        raise ValueError("Not all selected fields in table columns.")

    get_fields = op.itemgetter(
        *tuple(select_fields)
    )

    return get_fields


@fnt.lru_cache(4)
def decode_snv_allele_code(repr):

    if repr == "1:0:0:0":
        alt = "A"
    elif repr == "0:1:0:0":
        alt = "C"
    elif repr == "0:0:1:0":
        alt = "G"
    elif repr == "0:0:0:1":
        alt = "T"
    else:
        raise ValueError(f"Unknown SNV ALT repr: {repr}")
    return alt


def check_proper_merge_snv(stat_counter, get_fields, row):

    ref_pos = "\t".join([row["chrom"], row["start"], row["end"]])

    call_ids, samples, callsets, ref_alleles, alt_alleles = get_fields(row)
    call_ids = call_ids.split("|")
    samples = samples.split("|")
    callsets = callsets.split("|")
    alt_alleles = alt_alleles.split("|")
    alt_alleles = [decode_snv_allele_code(a) for a in alt_alleles]
    ref_alleles = ref_alleles.split("|")
    ref_alleles = [decode_snv_allele_code(r) for r in ref_alleles]

    out_info_singles = ""
    out_info_multis = ""

    if len(call_ids) > 1:
        # is multiple
        group_size = len(call_ids)
        group_id = hl.md5("".join(sorted(call_ids)).encode("utf-8")).hexdigest()
        for call_id, sample, callset, ref, alt in zip(call_ids, samples, callsets, ref_alleles, alt_alleles):
            stat_counter[(sample, callset, "SNV", "multiple")] += 1
            out_info_multis += (
                f"{ref_pos}\t{call_id}\tSNV\t1\t"
                f"{sample}\t{callset}\t"
                f"{group_id}\t{group_size}\t"
                f"{ref}\t{alt}\n"
            )
    else:
        # is singleton
        stat_counter[(samples[0], callsets[0], "SNV", "singleton")] += 1
        out_info_singles = (
            f"{ref_pos}\t{call_ids[0]}\t"
            f"SNV\t1\t{samples[0]}\t{callsets[0]}\t"
            f"{ref_alleles[0]}\t{alt_alleles[0]}\n"
        )
    return out_info_singles, out_info_multis


def check_proper_merge_indel(stat_counter, get_fields, row):

    ref_pos = "\t".join([row["chrom"], row["start"], row["end"]])

    call_ids, samples, callsets, var_types, var_lengths, ref_alleles, alt_alleles = get_fields(row)
    call_ids = call_ids.split("|")
    samples = samples.split("|")
    callsets = callsets.split("|")
    var_types = var_types.split("|")
    var_lengths = var_lengths.split("|")
    ref_alleles = ref_alleles.split("|")
    ref_ids = col.Counter(ref_alleles)
    ref_ids = dict(
        (ref_allele, allele_id) for allele_id, (ref_allele, allele_count)
        in enumerate(ref_ids.most_common(), start=1)
    )
    # for SYM calls, set this to zero
    ref_ids["0:0:0:0"] = 0

    alt_alleles = alt_alleles.split("|")
    alt_ids = col.Counter(alt_alleles)
    alt_ids = dict(
        (alt_allele, allele_id) for allele_id, (alt_allele, allele_count)
        in enumerate(alt_ids.most_common(), start=1)
    )
    alt_ids["0:0:0:0"] = 0

    out_info_singles = ""
    out_info_multis = ""

    if len(call_ids) > 1:
        # maybe multi-call ...
        # first group by variant type and length
        type_buckets = col.defaultdict(list)

        for var_type, var_length, call_id, sample, callset, ref, alt in zip(
                var_types, var_lengths, call_ids, samples, callsets, ref_alleles, alt_alleles
            ):

            type_buckets[(var_type, var_length)].append(
                (call_id, sample, callset, f"R{ref_ids[ref]}", f"A{alt_ids[alt]}")
            )

        for (var_type, var_length), calls in type_buckets.items():
            if len(calls) == 1:
                call_id, sample, callset, ref_repr, alt_repr = calls[0]
                stat_counter[(sample, callset, var_type, "singleton")] += 1
                out_info_singles = (
                    f"{ref_pos}\t{call_id}\t"
                    f"{var_type}\t{var_length}\t"
                    f"{sample}\t{callset}\t"
                    f"{ref_repr}\t{alt_repr}\n"
                )
            else:
                call_ids = [t[0] for t in calls]
                group_id = hl.md5("".join(sorted(call_ids)).encode("utf-8")).hexdigest()
                group_size = len(call_ids)

                samples = [t[1] for t in calls]
                callsets = [t[2] for t in calls]
                ref_reprs = [t[3] for t in calls]
                alt_reprs = [t[4] for t in calls]
                for call_id, sample, callset, ref_repr, alt_repr in zip(call_ids, samples, callsets, ref_reprs, alt_reprs):
                    stat_counter[(sample, callset, var_type, "multiple")] += 1
                    out_info_multis += (
                        f"{ref_pos}\t{call_id}\t"
                        f"{var_type}\t{var_length}\t"
                        f"{sample}\t{callset}\t"
                        f"{group_id}\t{group_size}\t"
                        f"{ref_repr}\t{alt_repr}\n"
                    )
    else:
        # is singleton
        stat_counter[(samples[0], callsets[0], var_types[0], "singleton")] += 1
        ref_repr = f"R{ref_ids[ref_alleles[0]]}"
        alt_repr = f"A{alt_ids[alt_alleles[0]]}"

        out_info_singles = (
            f"{ref_pos}\t{call_ids[0]}\t"
            f"{var_types[0]}\t{var_lengths[0]}\t"
            f"{samples[0]}\t{callsets[0]}\t"
            f"{ref_repr}\t{alt_repr}\n"
        )
    return out_info_singles, out_info_multis


def get_row_processor(variant_group):

    row_func = None
    if variant_group == "SNV":
        row_func = check_proper_merge_snv
    elif variant_group == "INDEL":
        row_func = check_proper_merge_indel
    else:
        raise NotImplementedError(f"No row processor for: {variant_group}")
    return row_func


def main():

    args = parse_command_line()

    # create folders for output
    args.singletons.parent.mkdir(exist_ok=True, parents=True)
    args.multicalls.parent.mkdir(exist_ok=True, parents=True)

    stat_counter = col.Counter()

    buffer_singles = io.StringIO()
    buffer_multis = io.StringIO()
    buffer_size = 0
    with ctl.ExitStack() as exs:
        merged_table = exs.enter_context(xopen.xopen(args.merged_table, "r"))

        column_names = check_or_get_column_header(merged_table, args.variant_group)
        get_columns = field_getter(args.variant_group, column_names)
        row_processor = get_row_processor(args.variant_group)
        row_processor = fnt.partial(row_processor, *(stat_counter, get_columns))

        singletons = exs.enter_context(xopen.xopen(args.singletons, "w"))
        singletons.write("\t".join(output_header("singletons")) + "\n")

        multicalls = exs.enter_context(xopen.xopen(args.multicalls, "w"))
        multicalls.write("\t".join(output_header("multicalls")) + "\n")

        reader = csv.DictReader(merged_table, fieldnames=column_names, delimiter="\t")
        for row in reader:
            out_singles, out_multis = row_processor(row)
            buffer_size += len(out_singles) + len(out_multis)
            buffer_singles.write(out_singles)
            buffer_multis.write(out_multis)

            if buffer_size > int(5e8):
                singletons.write(buffer_singles.getvalue())
                buffer_singles = io.StringIO()
                multicalls.write(buffer_multis.getvalue())
                buffer_multis = io.StringIO()

        if buffer_size > 0:
            singletons.write(buffer_singles.getvalue())
            multicalls.write(buffer_multis.getvalue())

    args.count_stats.parent.mkdir(exist_ok=True, parents=True)
    with open(args.count_stats, "w") as dump:
        _ = dump.write(f"sample\tcallset\tvar_type\tshared\tcount\n")
        for entity, count in stat_counter.most_common():
            row = "\t".join(entity) + f"\t{int(count)}\n"
            _ = dump.write(row)

    return 0


if __name__ == "__main__":
    main()
