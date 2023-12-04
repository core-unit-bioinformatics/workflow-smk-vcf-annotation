#!/usr/bin/env python3

import argparse as argp
import collections as col
import contextlib as ctl
import gzip
import hashlib as hl
import pathlib as pl

import numpy as np
import pandas as pd


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


def has_compatible_size(group_variant, test_variant):

    ref_length = group_variant.end - group_variant.start
    if ref_length < 2 and group_variant.size > 1:
        # no "length" in reference space, e.g., INS call
        # just compare sizes
        size_frac = test_variant.size / group_variant.size
        is_compatible = 0.75 < size_frac < 1.25
    else:
        # if length in reference space, check reciprocal overlap
        overlap = min(group_variant.end, test_variant.end) - max(group_variant.start, test_variant.start)
        ovl_group = overlap / group_variant.size
        ovl_test = overlap / test_variant.size
        is_compatible = ovl_group > 0.5 and ovl_test > 0.5
    return is_compatible


def check_grouping_criteria(sv_group, sv_to_check):

    same_type = sv_to_check.sv_merge_type == sv_group[0].sv_merge_type
    if not same_type:
        return False

    close_start = all(
        sv.start - 100 <= sv_to_check.start <= sv.start + 100
        for sv in sv_group
    )
    close_end = all(
        sv.end - 100 <= sv_to_check.end <= sv.end + 100
        for sv in sv_group
    )
    if not (close_start and close_end):
        return False

    recip_overlap = all(
        has_compatible_size(sv, sv_to_check) for sv in sv_group
    )
    if not recip_overlap:
        return False

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

    groups = []

    for row in sv_infos.itertuples():
        if not groups:
            groups.append([row])
            continue
        for num, sv_group in enumerate(groups, start=0):
            is_compatible = check_grouping_criteria(sv_group, row)
            if is_compatible:
                groups[num].append(row)
                break

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

    groups = dict()
    group_sizes = col.Counter()

    with open_merge(args.sv_merge, read_merge) as merge_table:
        for line in merge_table:
            chrom, start, end, concat_ids = line.strip().split()
            concat_ids = concat_ids.split("|")
            if len(concat_ids) == 1:
                continue
            else:
                # NB: identical calls can have the same name
                # (= same hash), which necessitates making
                # the list of names unique to avoid selecting
                # calls several times
                concat_ids = list(set(concat_ids))
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
