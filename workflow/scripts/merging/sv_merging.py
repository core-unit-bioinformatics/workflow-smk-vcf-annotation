#!/usr/bin/env python3

import argparse as argp
import collections as col
import contextlib as ctl
import gzip
import hashlib as hl
import pathlib as pl
import sys

import numpy as np
import pandas as pd


DEBUG_MODE = False


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--sv-table",
        "-st",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="sv_table",
        help="Path to TSV table with complete SV information",
        required=True
    )

    parser.add_argument(
        "--sv-merge",
        "-sm",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="sv_merge",
        help="Path to TSV with positionally (bedtools) merged SVs.",
        required=True
    )

    parser.add_argument(
        "--debug-ids",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="debug_ids",
        help="Path to text file listing variant IDs for a debug run. Default: None",
        default=None
    )

    parser.add_argument(
        "--singletons",
        "-s", "-sng",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="singletons",
        help="Path to output file for singleton SVs.",
        required=True
    )

    parser.add_argument(
        "--multiples",
        "-m", "-mul",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="multiples",
        help="Path to output file for shared SVs.",
        required=True
    )

    parser.add_argument(
        "--count-stats",
        "-stats", "-c",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="count_stats",
        help="Path to output file for count statistics.",
        required=True
    )

    parser.add_argument(
        "--desc-stats",
        "-desc", "-d",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="desc_stats",
        help="Path to output file for SV characteristics / descriptive stats.",
        required=True
    )

    args = parser.parse_args()
    return args


def get_sv_table_header():
    """Ugly hack ... read the header information
    from the "Variant" class definition in the
    tab_size_convert.py script - horrible ...
    """
    my_loc = pl.Path(__file__).resolve(strict=True)
    script_loc = my_loc.parent.parent.joinpath("table_convert", "tab_size_dist.py")
    assert script_loc.is_file(), f"Cannot extract table header from script: {script_loc}"

    table_header = []
    in_class_def = False
    fetch_info = False
    with open(script_loc, "r") as script_file:
        for line in script_file:
            if line.startswith("class Variant:"):
                in_class_def = True
                continue
            if "__slots__ = (" in line.strip() and in_class_def:
                fetch_info = True
                continue
            if ")" in line and in_class_def:
                break
            if fetch_info:
                column_name = line.strip().split()[0]
                column_name = column_name.strip("\",").strip()
                assert column_name.isidentifier(), f"Invalid column name: {column_name}"
                table_header.append(column_name)

    return table_header


def determine_file_type(file_path):

    open_func = open
    write_mode = "w"
    read_mode = "r"
    if file_path.suffix == ".gz":
        open_func = gzip.open
        write_mode = "wt"
        read_mode = "rt"
    return open_func, read_mode, write_mode


def has_compatible_spread(group_variant, test_variant, lenient_end_check):

    global DEBUG_MODE

    ref_length = group_variant.end - group_variant.start
    if (ref_length < 2 and group_variant.size > 1) or lenient_end_check:
        # no "length" in reference space, e.g., INS call
        # just compare sizes
        size_frac = test_variant.size / group_variant.size
        if DEBUG_MODE:
            sys.stderr.write(f"\nHCS - case 1 / size_frac: {size_frac}\n")
        is_compatible = 0.75 < size_frac < 1.25
    elif (ref_length < 2 and group_variant.size < 2):
        assert group_variant.sv_merge_type in ["BND"]
        if DEBUG_MODE:
            sys.stderr.write(f"\nHCS - case 2 / BND\n")
        # BNDs are compatible iff
        # chrom2 is identical
        # chrom2_pos is close
        # join operation is identical
        # === important to realize that lenient_end_check is only
        # === triggered for INS/DUP, those do not end up here
        chrom2_is_identical = group_variant.chrom2 == test_variant.chrom2
        chrom2_pos_is_close = test_variant.chrom2_pos - 100 <= group_variant.chrom2_pos <= test_variant.chrom2_pos + 100
        join_op_is_identical = group_variant.bnd_join_vcf == test_variant.bnd_join_vcf
        is_compatible = chrom2_is_identical & chrom2_pos_is_close & join_op_is_identical
    else:
        # if length in reference space, check reciprocal overlap
        overlap = min(group_variant.end, test_variant.end) - max(group_variant.start, test_variant.start)
        if DEBUG_MODE:
            sys.stderr.write(f"\nHCS - case 3 / ref_length {ref_length} || overlap {overlap}\n")
        ovl_group = overlap / group_variant.size
        ovl_test = overlap / test_variant.size
        is_compatible = ovl_group > 0.5 and ovl_test > 0.5
    return is_compatible


def check_grouping_criteria(sv_group, sv_to_check):

    global DEBUG_MODE
    DEBUG_RETURN = True  # records what should have been returned

    same_type = sv_to_check.sv_merge_type == sv_group[0].sv_merge_type
    if DEBUG_MODE:
        DEBUG_RETURN &= same_type
    elif not same_type:
        return False
    else:
        pass

    close_start = all(
        sv.start - 100 <= sv_to_check.start <= sv.start + 100
        for sv in sv_group
    )
    close_end = all(
        sv.end - 100 <= sv_to_check.end <= sv.end + 100
        for sv in sv_group
    )

    # Annoying special case: a caller like pbsv produces DUPs
    # with two breakpoints in reference space whereas other callers
    # may just call an INS with a single breakpoint, which means
    # the end breakpoint will not generally fall within the distance
    # even if the two calls are likely the same event.
    # If and only if the SV types belong to these two types,
    # the "close_end" criterion can be ignored/relaxed.
    # TODO: turn this into a config parameter?
    lenient_end_check = False
    if not (close_start and close_end):
        # the actual variant type of a group can be mixed for the INS/DUP
        # case, hence the line below selects for the majority type
        actual_group_types = set(sv.vartype for sv in sv_group)
        group_is_homogen = len(actual_group_types) == 1
        group_type = None
        if group_is_homogen:
            group_type = actual_group_types.pop()
        actual_check_type = sv_to_check.vartype

        # Important for debugging: after merging the first DUP call into a group of INS,
        # the all() checks for the coordinates above will most likely fail. Need to consider
        # a few scenarios in the following:

        # if the group is homogeneous and the SV type to check is different, only
        # INS/DUP combinations can be checked leniently
        if group_is_homogen and actual_check_type != group_type:
            lenient_end_check = group_type == "DUP" and actual_check_type == "INS"
            lenient_end_check |= group_type == "INS" and actual_check_type == "DUP"
        elif group_is_homogen:
            # implies actual_check_type is identical to group type, hence the
            # end coordinates should match; maybe, e.g., different ALT here
            lenient_end_check = False
        else:
            # group is not homogeneous, e.g., DUP already merged into INS
            # if so, lenient check if everything is DUP/INS
            lenient_end_check = actual_check_type in ["DUP", "INS"]
            lenient_end_check &= all(grp in ["DUP", "INS"] for grp in actual_group_types)

        if close_start and lenient_end_check:
            # could merge, RO-test may or may not fail
            pass
        elif DEBUG_MODE:
            DEBUG_RETURN &= (close_start and lenient_end_check)
        else:
            return False

    recip_overlap = all(
        has_compatible_spread(sv, sv_to_check, lenient_end_check) for sv in sv_group
    )
    if DEBUG_MODE:
        DEBUG_RETURN &= recip_overlap
    elif not recip_overlap:
        return False
    else:
        pass

    if DEBUG_MODE:
        DEBUG_CHECKS = [
            ("same_type", same_type),
            ("lenient_end_check", lenient_end_check),
            ("recip_overlap", recip_overlap),
            ("close_start", close_start),
            ("close_end", close_end)
        ]
        try:
            _ = group_is_homogen
        except UnboundLocalError:
            group_is_homogen = "untested"
        DEBUG_CHECKS.append(("group_is_homogen", group_is_homogen))
        sys.stderr.write("\nDEBUG print ---\n")
        sys.stderr.write(" || ".join([f"{name}: {value}" for name, value in DEBUG_CHECKS]) + "\n")
        sys.stderr.write(" --- end\n")
        return DEBUG_RETURN

    return True



def check_sv_grouping(sv_infos):
    """Check if overlapping SV calls
    can be grouped (likely represent the same event).

    NB: this function assumes that only calls
    on the same reference chromosome are being
    compared, so that is not checked.

    Args:
        sv_infos (_type_): _description_
    """
    global DEBUG_MODE
    groups = []

    for row in sv_infos.itertuples():
        if not groups:
            groups.append([row])
            continue
        is_first_member = True
        for num, sv_group in enumerate(groups, start=0):
            is_compatible = check_grouping_criteria(sv_group, row)
            if is_compatible:
                is_first_member = False
                groups[num].append(row)
                break

        if is_first_member:
            groups.append([row])

    call_to_group = dict()
    group_to_size = dict()
    for group in groups:
        group_size = len(group)
        if group_size == 1:
            # could happen for first row added to groups
            continue
        call_ids = sorted(g.Index for g in group)
        group_id = hl.md5("".join(call_ids).encode("utf-8")).hexdigest()
        call_to_group.update(
            dict((call_id, group_id) for call_id in call_ids)
        )
        group_to_size[group_id] = group_size

    if DEBUG_MODE:
        sys.stderr.write("\nDEBUG print ---\n")
        sys.stderr.write(" || ".join([f"GRP id {group} -- size {size}" for group, size in group_to_size.items()]) + "\n")
        sys.stderr.write(" --- end\n")

    return call_to_group, group_to_size


def dump_counting_statistics(out_file, count_stats):

    out_file.parent.mkdir(exist_ok=True, parents=True)
    with open(out_file, "w") as dump:
        _ = dump.write(f"sample\tcallset\tvar_type\tshared\tcount\n")
        for row in count_stats.itertuples():
            sample, callset, vartype = row.Index
            dump.write(f"{sample}\t{callset}\t{vartype}\tsingleton\t{row.singleton}\n")
            dump.write(f"{sample}\t{callset}\t{vartype}\tmultiple\t{row.multiple}\n")
    return


def dump_variant_calls(out_file, subset, sv_infos):

    common_header = [
        "chrom", "start", "end", "name", "vartype", "size",
        "sample", "callset"
    ]

    if subset == "singleton":
        header = common_header + [
             "ref_allele_repr", "alt_allele_repr"
        ]
        selected_calls = sv_infos.loc[sv_infos["group_id"] == "singleton", :]
    elif subset == "multiple":
        header = common_header + [
            "group_id", "group_size", "ref_allele_repr", "alt_allele_repr"
        ]
        selected_calls = sv_infos.loc[sv_infos["group_id"] != "singleton", :]
    else:
        raise

    out_file.parent.mkdir(exist_ok=True, parents=True)
    selected_calls[header].to_csv(out_file, sep="\t", header=True, index=False)

    return


def extract_sv_statistics(sv_infos):

    sv_infos["is_shared"] = 0
    sv_infos.loc[sv_infos["group_id"] != "singleton", "is_shared"] = 1

    select_qual_stats = [
        "size", "read_depth", "call_quality", "genotype_quality",
        "ref_allele_freq", "alt_allele_freq"
    ]

    norm_stat_name = {
        "count": "count",
        "mean": "mean",
        "std": "stddev",
        "min": "min",
        "max": "max"
    }
    norm_stat_name.update(
        dict((f"{n}%", f"pctile_{n}") for n in [5,10,25,50,75,90,95])
    )

    merge_stats = []
    for (vartype, is_shared), stats in sv_infos.groupby(["vartype", "is_shared"])[select_qual_stats]:
        sub = stats.replace(to_replace=-1, value=np.nan)
        desc_stats = sub.describe(percentiles=[0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95])
        category = "singleton" if is_shared < 1 else "multiple"
        new_index = pd.MultiIndex.from_tuples(
            [(vartype, category, norm_stat_name[c]) for c in desc_stats.index.values],
            names=["vartype", "occurrence", "statistic"]
        )
        desc_stats.set_index(new_index, inplace=True)
        merge_stats.append(desc_stats)

    merge_stats = pd.concat(merge_stats, axis=0, ignore_index=False)
    merge_stats = merge_stats.round(5)

    return merge_stats


def main():

    global DEBUG_MODE

    args = parse_command_line()
    sv_table_header = get_sv_table_header()

    # this is one table (per chromosome) for all samples/callsets,
    # so this structure can be huge in memory
    sv_infos = pd.read_csv(
        args.sv_table, sep="\t", comment="#",
        header=None, names=sv_table_header,
        index_col="name", low_memory=False
    )
    # empirically, callers may differ in the way they are labeling
    # INS/DUPs, so recode here for easier merge later on
    sv_merge_types = {
        "DEL": "DEL", "INS": "INS", "DUP": "INS",
        "BND": "BND", "INV": "INV"
    }
    assert all(uniq in sv_merge_types for uniq in sv_infos["vartype"].unique()), \
        f"{sv_infos['vartype'].unique()}"

    sv_infos["sv_merge_type"] = sv_infos["vartype"].replace(sv_merge_types)

    open_merge, read_merge, _ = determine_file_type(args.sv_merge)

    if args.debug_ids is not None:

        DEBUG_MODE = True
        with open(args.debug_ids, "r") as listing:
            DEBUG_IDS = listing.read().strip().split()
        DEBUG_PROCESSED = False

    groups = dict()
    group_sizes = col.Counter()

    with open_merge(args.sv_merge, read_merge) as merge_table:
        for line in merge_table:
            if DEBUG_MODE and DEBUG_PROCESSED:
                break
            chrom, start, end, concat_ids = line.strip().split()
            concat_ids = concat_ids.split("|")
            if DEBUG_MODE:
                if any(debug_id in concat_ids for debug_id in DEBUG_IDS):
                    concat_ids = DEBUG_IDS
                    DEBUG_PROCESSED = True
                else:
                    continue
            if len(concat_ids) == 1:
                continue
            else:
                # NB: identical calls can have the same name
                # (= same hash), which necessitates making
                # the list of names unique to avoid selecting
                # calls several times
                concat_ids = sorted(set(concat_ids))
                group_sv_infos = sv_infos.loc[concat_ids, :]
                grouped_calls, sizes = check_sv_grouping(group_sv_infos)
                groups.update(grouped_calls)
                group_sizes.update(sizes)

    sv_infos["group_id"] = sv_infos.index.map(lambda x: groups.get(x, "singleton"))
    sv_infos["group_size"] = sv_infos["group_id"].apply(lambda x: group_sizes[x])
    sv_infos = sv_infos.reset_index(drop=False, inplace=False)

    count_stats = sv_infos.groupby(["sample", "callset", "vartype"])["group_id"].agg(
        singleton=lambda x: (x == "singleton").sum(),
        multiple=lambda x: (x != "singleton").sum()
    )
    dump_counting_statistics(args.count_stats, count_stats)

    desc_stats = extract_sv_statistics(sv_infos)
    args.desc_stats.parent.mkdir(exist_ok=True, parents=True)
    desc_stats.to_csv(
        args.desc_stats, sep="\t",
        header=True, index=True, na_rep="nan"
    )

    dump_variant_calls(args.singletons, "singleton", sv_infos)
    dump_variant_calls(args.multiples, "multiple", sv_infos)

    return 0


if __name__ == "__main__":
    main()
