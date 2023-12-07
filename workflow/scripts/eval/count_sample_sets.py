#!/usr/bin/env python3

import argparse as argp
import collections as col
import itertools as itt
import pathlib as pl

import pandas as pd
import xopen as xo


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--grouped-calls", "-g",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="grouped_calls",
        help="Path to TSV table with grouped variant calls and sample/callset indicator columns."
    )

    parser.add_argument(
        "--out-subsets", "-sub",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_subsets",
        help="Path to TSV output table listing all sample/callset subsets sharing a call."
    )

    parser.add_argument(
        "--out-subset-call-map", "-map",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_mapping",
        help="Path to TSV output table mapping sample / callset groups to call group ID."
    )

    args = parser.parse_args()

    return args


def determine_sample_callset_columns(table_columns):

    info_columns = [
        "chrom", "start", "end", "size",
        "vartype", "group_size", "distinct_samples"
    ]

    if not all(c in table_columns for c in info_columns):
        raise ValueError(
            "Not all call/group information column names in table header. "
            "Did the input format change?\n"
            f"Info columns: {info_columns}\n"
            f"Table header: {table_columns}\n"
        )

    sample_callset_columns = [c for c in table_columns if c not in info_columns]

    return info_columns, sample_callset_columns


def count_sample_callset_subsets(call_table, indicator_columns):

    sample_subsets = col.Counter()
    callset_subsets = col.Counter()

    sets_to_group = col.defaultdict(set)

    for idx, row in call_table.iterrows():
        # 1 - drop all sample/callsets that are n/a
        sample_callsets = row[indicator_columns].dropna(inplace=False).index
        # 2 - split into samples and callsets
        # === TODO: split by "_" is a hard-coded delimiter following
        # === the script merge_by_group ; not optimal ...
        # Note: for samples, counting duplicates makes no sense because
        # these would just be calls from different callsets. For callsets,
        # however, it does make sense because the same callset identifier
        # belongs to different samples. Hence, counting these duplicates
        # would potentially reveal calls that are only made, e.g., by one
        # caller but none of the others (or combination of aligner and caller
        # and so on)
        samples = sorted(set([entry.split("_", 1)[0] for entry in sample_callsets]))
        callsets = sorted([entry.split("_", 1)[-1] for entry in sample_callsets])

        sample_subsets[tuple(samples)] += 1
        sets_to_group[tuple(samples)].add(idx)
        if len(samples) > 2:
            for s1, s2 in itt.combinations(samples, r=2):
                sample_subsets[(s1, s2)] += 1
                sets_to_group[(s1, s2)].add(idx)

        callset_subsets[tuple(callsets)] += 1
        sets_to_group[tuple(callsets)].add(idx)
        if len(callsets) > 2:
            for c1, c2 in itt.combinations(callsets, r=2):
                callset_subsets[(c1, c2)] += 1
                sets_to_group[(c1, c2)].add(idx)

    return sample_subsets, callset_subsets, sets_to_group


def main():

    args = parse_command_line()

    call_table = pd.read_csv(
        args.grouped_calls, sep="\t",
        header=0, comment="#",
        low_memory=False
    )
    call_table.set_index("group_id", inplace=True)

    info_cols, sample_callset_cols = determine_sample_callset_columns(call_table.columns)

    sample_subsets, callset_subsets, sets_to_callgroup = count_sample_callset_subsets(
        call_table, sample_callset_cols
    )

    args.out_subsets.parent.mkdir(exist_ok=True, parents=True)
    with xo.xopen(args.out_subsets, "w") as dump:
        _ = dump.write("category\tcardinality\tsubset\tabundance\n")
        for comb, count in sample_subsets.most_common():
            _ = dump.write(f"samples\t{len(comb)}\t{','.join(comb)}\t{count}\n")
        for comb, count in callset_subsets.most_common():
            _ = dump.write(f"callsets\t{len(comb)}\t{','.join(comb)}\t{count}\n")

    args.out_mapping.parent.mkdir(exist_ok=True, parents=True)
    with xo.xopen(args.out_mapping, "w") as dump:
        _ = dump.write("subset\tcall_group_ids\n")
        for subset, call_ids in sets_to_callgroup.items():
            _ = dump.write(f"{','.join(subset)}\t{','.join(sorted(call_ids))}\n")

    return 0


if __name__ == "__main__":
    main()
