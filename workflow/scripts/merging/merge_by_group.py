#!/usr/bin/env python3

import argparse as argp
import pathlib as pl


import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--input-table",
        "-i", "-tab",
        nargs="+",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="input_table",
        help="Path to TSV input tables with grouped variant calls."
    )

    parser.add_argument(
        "--output-table",
        "-o", "-tsv",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="output_table",
        help=(
            "Path to TSV output table. If several inputs are provided, "
            "they will be concatenated in the output."
        )
    )

    args = parser.parse_args()

    return args


def make_by_group_table(tab_calls):

    tab_calls["sample_callset"] = tab_calls["sample"] + "_" + tab_calls["callset"]

    sample_callset_columns = set(tab_calls["sample_callset"].unique())

    by_group = tab_calls.groupby("group_id").agg(
        {
            "chrom": pd.Series.mode,
            "start": "min",
            "end": "max",
            "size": lambda s: int(s.median()),
            "vartype": lambda s: s.mode()[0],  # in case of ties, select first
            "group_size": "max"
        }
    )
    by_group["distinct_samples"] = 0

    # indicator = indicate which calls/groups
    # appear in which callset
    indicator = pd.DataFrame(
        [], dtype=str, index=by_group.index,
        columns=tab_calls["sample_callset"].unique()
    )
    indicator.fillna("n/a", inplace=True)

    for group_id, sample_calls in tab_calls.groupby("group_id"):
        # TODO - for some 'sample_calls' objects, ".pivot" was raising
        # because of the duplicated entries in index (= group id). Unclear
        # why it worked for a subset w/o issues. Hence, this slow and
        # manual solution to building the call-to-sample/callset indicator
        by_group.loc[group_id, "distinct_samples"] = sample_calls["sample"].nunique()
        for row in sample_calls.itertuples():
            indicator.loc[row.group_id, row.sample_callset] = row.name

    by_group = by_group.join(indicator, how="outer")
    group_ids = by_group.index.values
    by_group.reset_index(drop=True, inplace=True)
    by_group.insert(3, column="group_id", value=group_ids)

    return by_group, sample_callset_columns


def main():

    args = parse_command_line()

    concat_calls = []
    indicator_columns = set()
    for table_file in args.input_table:
        tab_calls = pd.read_csv(table_file, sep="\t", header=0, comment="#")
        if tab_calls.empty:
            continue
        grouped_calls, callset_columns = make_by_group_table(tab_calls)
        indicator_columns = indicator_columns.union(callset_columns)
        concat_calls.append(grouped_calls)

    indicator_columns = sorted(indicator_columns)
    concat_calls = pd.concat(concat_calls, axis=0, ignore_index=False)
    concat_calls[indicator_columns] = concat_calls[indicator_columns].fillna("n/a", inplace=False)
    concat_calls.sort_values(["chrom", "start", "end"], inplace=True)

    args.output_table.parent.mkdir(exist_ok=True, parents=True)
    concat_calls.to_csv(args.output_table, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
