#!/usr/bin/env python3

import argparse as argp
import contextlib as ctl
import hashlib as hl
import pathlib as pl
import json
import sys

import pandas as pd
import xopen


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--input",
        "-i",
        type=str,
        default="stdin",
        dest="input",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="stdout",
        dest="output"
    )

    parser.add_argument(
        "--name-idx",
        "-n",
        type=lambda x: pl.Path(x).resolve(),
        default="vcf_name.idx.tsv.gz",
        dest="idx_tsv"
    )

    parser.add_argument(
        "--idx-info",
        "-j",
        type=lambda x: pl.Path(x).resolve(),
        default="vcf_name.idx.info.json",
        dest="idx_json"
    )

    parser.add_argument(
        "--tags",
        "-t",
        nargs="*",
        type=str,
        dest="tags"
    )

    args = parser.parse_args()

    return args


def main():

    args = parse_command_line()

    name_index = []
    record_num = 0
    input_filename = None
    with ctl.ExitStack() as stack:
        if args.input == "stdin":
            read_input = sys.stdin
        else:
            input_file = pl.Path(args.input)
            assert input_file.is_file()
            input_filename = input_file.name
            read_input = stack.enter_context(xopen.xopen(args.input))
        if args.output == "stdout":
            write_output = sys.stdout
        else:
            outfile = pl.Path(args.output).resolve()
            outfile.parent.mkdir(exist_ok=True, parents=True)
            write_output = stack.enter_context(xopen.xopen(outfile, "wt"))

        for ln, line in enumerate(read_input, 1):
            if line.startswith("#"):
                write_output.write(line)
                continue
            record_num += 1
            new_name = hl.md5(line.strip().encode("utf-8")).hexdigest()
            columns = line.strip().split()
            old_name = columns[2]
            columns[2] = new_name
            write_output.write("\t".join(columns) + "\n")
            name_index.append((ln, record_num, new_name, old_name))

        stack.close()

    name_index = pd.DataFrame.from_records(
        name_index, columns=["line_num", "record_num", "idx_name", "name"]
    )
    assert record_num == name_index.shape[0]
    args.idx_tsv.parent.mkdir(exist_ok=True, parents=True)
    name_index.to_csv(args.idx_tsv, sep="\t", header=True, index=False)

    idx_info = {
        "vcf_records": record_num,
    }
    if input_filename is not None:
        idx_info["filename"] = input_filename

    if args.tags:
        for tag in args.tags:
            key, value = tag.split(":", 1)
            idx_info[key] = value

    args.idx_json.parent.mkdir(exist_ok=True, parents=True)
    with open(args.idx_json, "w") as dump:
        _ = json.dump(idx_info, dump)

    return 0


if __name__ == "__main__":
    sys.exit(main())
