#!/usr/bin/env python3

import argparse as argp
import pathlib as pl

import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--positive", "--positive-tables", "-p",
        type=lambda x: pl.Path(x).resolve(strict=True),
        nargs="+",
        dest="positive_tables"
    )

    parser.add_argument(
        "--negative", "--negative-tables", "-n",
        type=lambda x: pl.Path(x).resolve(strict=True),
        nargs="*",
        dest="negative_tables"
    )

    parser.add_argument(
        "--distance-criterion", "-d",
        type=str,
        choices=["fixed", "relative"],
        default="relative",
        dest="distance_criterion"
    )

    parser.add_argument(
        "--distance-cutoff", "-c",
        type=float,
        dest="distance_cutoff",
        default=0.1
    )

    parser.add_argument(
        "--output", "-o",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="output",
        default=pl.Path("joined.tsv.gz")
    )

    parser.add_argument(
        "--size-limit", "-l",
        type=int,
        default=int(1e5),
        dest="size_limit"
    )

    args = parser.parse_args()

    return args


def read_table(filepath, size_limit):

    df = pd.read_csv(filepath, sep="\t", header=0)

    # we hard drop spurious calls that are very
    # large (by default: 100+ kbp)
    df = df.loc[df["size"] < size_limit, :].copy()

    dist_column = df.columns[-1]
    assert dist_column.startswith("distance_")
    annotation_name = dist_column.split("_")[-1]

    # fix: bedtools intersect may create a dist of -1
    # for non-matching regions --- drop those lines immediately
    df = df.loc[df[dist_column] >= 0, :].copy()

    # some columns that can always be dropped
    hard_drop = ["chrom", "start", "end", "score", "strand", "chrom"]
    drop_from_ann = [f"{column}_{annotation_name}" for column in hard_drop]

    # heuristic
    drop_from_ann = drop_from_ann + [
        c for c in df.columns if "_start_" in c or "_end_" in c
    ]

    df.drop(drop_from_ann, axis=1, inplace=True, errors="ignore")

    return df



def get_deselect_call_ids(negative_table, dist_criterion, dist_cutoff):

    dist_column = negative_table.columns[-1]

    if dist_criterion == "relative":
        # relative cutoff
        # DISTANCE < SIZE_OF_CALL * DIST_CUTOFF
        # Example default values:
        # 1000 * 0.1 => 100
        # leads to:
        # DISTANCE < 100 => discard call
        # if there is an entry in a negative annotation table
        # closer than 100 bp, the call is discarded
        # (the call group ID is included in the discard set)
        assert dist_cutoff < 1
        deselect = negative_table.loc[
            negative_table[dist_column] < negative_table["size"] * dist_cutoff,
            "group_id"
        ].values
    else:
        assert dist_cutoff > 1
        deselect = negative_table.loc[
            negative_table[dist_column] < dist_cutoff,
            "group_id"
        ].values

    return set(deselect)


def build_negative_set(negative_tables, dist_criterion, dist_cutoff, size_limit):

    deselect_calls = set()

    for negative_table in negative_tables:
        df_neg = read_table(negative_table, size_limit)
        delesect_from_table = get_deselect_call_ids(df_neg, dist_criterion, dist_cutoff)
        deselect_calls = deselect_calls.union(delesect_from_table)

    return deselect_calls


def main():

    args = parse_command_line()

    deselect_calls = build_negative_set(
        args.negative_tables, args.distance_criterion, args.distance_cutoff, args.size_limit
    )

    joined = None
    for annotation_table in sorted(args.positive_tables):
        ann_data = read_table(annotation_table, args.size_limit)
        ann_data = ann_data.loc[
            ~ann_data["group_id"].isin(deselect_calls),
            :
        ].copy()

        if ann_data.empty:
            continue

        dist_column = ann_data.columns[-1]
        if args.distance_criterion == "relative":
            # inverse to negative / discard
            # keep calls relatively close to an annotated region
            # such as a gene or an enhancer
            assert args.distance_cutoff < 1
            ann_data = ann_data.loc[
                ann_data[dist_column] < (ann_data["size"] * args.distance_cutoff),
                :
            ].copy()
        else:
            assert args.distance_cutoff > 1
            ann_data = ann_data.loc[
                ann_data[dist_column] < args.distance_cutoff,
                :
            ].copy()

        if joined is None:
            joined = ann_data
            continue
        common_columns = set(ann_data.columns).intersection(set(joined.columns))
        joined = joined.merge(
            ann_data, how="outer", left_index=False, right_index=False,
            on=list(common_columns)
        )

    args.output.parent.mkdir(exist_ok=True, parents=True)
    joined.sort_values(["#chrom", "start", "end"], inplace=True)

    # drop all columns that are empty
    empty_columns = joined.columns[pd.isnull(joined).all(axis=0).values]
    if len(empty_columns) > 0:
        joined.drop(empty_columns, axis=1, inplace=True)

    joined.to_csv(args.output, sep="\t", header=True, index=False, na_rep='n/a')

    return 0


if __name__ == "__main__":
    main()
