#!/usr/bin/env python3

import argparse as argp
import pathlib as pl

import numpy as np
import pandas as pd
import xopen


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        "--read-header", "-rhd",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="read_header",
        default=None,
        help=(
            "Read table header from this file (first line). "
            "If set to None (the default), the input table "
            "must have a valid header."
        )
    )

    parser.add_argument(
        "--call-table", "-calls", "-t",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="call_table",
        required=True,
        help="The table call table of all variant calls to subset."
    )

    parser.add_argument(
        "--singletons", "-sng", "-s",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="singletons",
        required=True,
        help="The table containing singleton calls."
    )

    parser.add_argument(
        "--stats-bracket", "-b",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="stats_bracket",
        required=True,
        help="The statistics table defining which singletons to select."
    )

    parser.add_argument(
        "--out-table", "-ot",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_table",
        required=True,
        help="Output table to dump selected singleton calls."
    )

    parser.add_argument(
        "--out-bed", "-ob",
        type=lambda x: pl.Path(x).resolve(strict=False),
        default=None,
        dest="out_bed",
        help="If set, output reduced table in a BED-like format. Default: None"
    )

    parser.add_argument(
        "--empty-output", "-empty",
        action="store_true",
        default=False,
        dest="empty_output",
        help=(
            "If no calls are selected, nevertheless "
            "create a header-only output table. "
            "Default: False"
        )
    )

    args = parser.parse_args()
    return args


def read_header(file_path):

    with xopen.xopen(file_path) as table:
        first_line = table.readline()
        if "," in first_line:
            header = first_line.strip().split(",")
        else:
            header = first_line.strip().split()
    assert len(header) > 1
    return header


def main():

    args = parse_command_line()

    singletons = set(
        pd.read_csv(
            args.singletons, sep="\t", header=0, usecols=["name"]
        )["name"].values
    )

    header = read_header(args.read_header)
    call_table = pd.read_csv(
        args.call_table, sep="\t", header=None,
        names=header, low_memory=False
    )
    call_table = call_table.loc[call_table["name"].isin(singletons), :].copy()

    stats_bracket = pd.read_csv(args.stats_bracket, sep="\t", header=0)
    stats_bracket.drop("occurrence", axis=1, inplace=True)

    vartypes = set(stats_bracket["vartype"].values).intersection(
        set(call_table["vartype"].values)
    )

    shared_stats = [
        c for c in stats_bracket.columns if c in call_table.columns
        and c != "vartype"
    ]
    assert len(shared_stats) > 1

    selector = np.ones(call_table.shape[0], dtype=bool)
    selected_calls = set()
    for vartype in vartypes:
        selector &= call_table["vartype"] == vartype
        var_limits = stats_bracket.loc[stats_bracket["vartype"] == vartype, :]
        assert var_limits.shape[0] == 2
        for s in shared_stats:
            selector &= call_table[s] >= var_limits[s].iloc[0]
            selector &= call_table[s] <= var_limits[s].iloc[1]
        if selector.any():
            selected_calls = selected_calls.union(
                set(call_table.loc[selector, "name"].values)
            )
        selector[:] = True

    if selected_calls:
        subset = call_table.loc[call_table["name"].isin(selected_calls), :].copy()

        args.out_table.parent.mkdir(exist_ok=True, parents=True)
        subset.to_csv(args.out_table, sep="\t", header=True, index=False)

        if args.out_bed is not None:
            args.out_bed.parent.mkdir(exist_ok=True, parents=True)
            with xopen.xopen(args.out_bed, "w") as bed_like:
                _ = bed_like.write("#")
                subset[
                    ["chrom", "start", "end", "name", "size", "vartype"]
                ].to_csv(bed_like, sep="\t", header=True, index=False)

    elif args.empty_output:
        args.out_table.parent.mkdir(exist_ok=True, parents=True)
        with xopen.xopen(args.out_table, "w") as dump:
            _ = dump.write("\t".join(header) + "\n")

        if args.out_bed is not None:
            args.out_bed.parent.mkdir(exist_ok=True, parents=True)
            with xopen.xopen(args.out_bed, "w") as bed_like:
                _ = bed_like.write("#")
                _ = bed_like.write("\t".join(header) + "\n")
    else:
        pass

    return 0


if __name__ == "__main__":
    main()
