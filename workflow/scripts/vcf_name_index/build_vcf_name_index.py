#!/usr/bin/env python3

import argparse as argp
import collections as col
import contextlib as ctl
import hashlib as hl
import itertools as itt
import pathlib as pl
import re
import json
import sys

import pandas as pd
import xopen


# global DEBUG switch
# currently only used for
# BDN partner renaming
DEBUG = False

VCF_INFO_FIELD_INDEX = 7


BREAKENDS = col.namedtuple("BREAKENDS", ["chrom1", "pos1", "chrom2", "pos2", "ln"])


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

    parser.add_argument(
        "--debug", "-dbg",
        action="store_true",
        default=False
    )

    args = parser.parse_args()

    if args.debug:
        global DEBUG
        DEBUG = True

    return args


def debug_msg(msg):
    full_msg = f"\nDEBUG:\n{msg}\n"
    sys.stderr.write(full_msg)
    return


def parse_vcf_header(header_lines):

    is_sv_callset = False
    has_mateid_field = False
    known_chroms = set()
    for line in header_lines:
        if line.startswith("##contig=<ID="):
            chrom = line.split(",")[0].split("=")[-1].strip()
            known_chroms.add(chrom)
            continue
        if not line.startswith("##INFO"):
            continue
        if "ID=SVTYPE" in line:
            is_sv_callset = True
            if DEBUG:
                debug_msg("Detected SV callset")
        if "ID=MATEID" in line:
            has_mateid_field = True
            if DEBUG:
                debug_msg("SV callset has MATEID info")
    # by reverse-sorting chromosome names by length, we get the longest match first
    match_chrom = re.compile("(" + "|".join(sorted(known_chroms, key=lambda x: len(x), reverse=True)) + "):[0-9]+")
    return is_sv_callset, has_mateid_field, match_chrom


def rename_exactly_by_mateid(bnd_lines):
    """For SV callers that report a MATEID such as pbsv
    """

    collect_partners = col.defaultdict(set)
    group_members = col.defaultdict(set)
    member_to_group_id = dict()
    group_counter = 0

    for (ln, rn, bnd_line) in bnd_lines:
        call_id = bnd_line[2]
        try:
            group_id = member_to_group_id[call_id]
        except KeyError:
            group_counter += 1
            group_id = group_counter
            group_members[group_id].add(call_id)
            member_to_group_id[call_id] = group_id

        mate_found = False
        for entry in bnd_line[VCF_INFO_FIELD_INDEX].split(";"):
            if entry.startswith("MATEID="):
                mate_ids = entry.replace("MATEID=", "").strip()
                for mate_id in mate_ids.split(","):
                    collect_partners[group_id].add((ln, rn, call_id, mate_id, bnd_line[0], tuple(bnd_line)))
                    group_members[group_id].add(mate_id)
                    if mate_id in member_to_group_id:
                        assert member_to_group_id[mate_id] == group_id
                    else:
                        member_to_group_id[mate_id] = group_id
                    mate_found = True
                    if DEBUG:
                        debug_msg(f"Paired: {call_id} / {mate_id}")
                break

        if not mate_found:
            if DEBUG:
                debug_msg(f"Singleton BND --- {bnd_line}")
            collect_partners[group_id].add((ln, rn, call_id, None, bnd_line[0], tuple(bnd_line)))
            # should be singleton BND

    return collect_partners, group_members


def rename_fuzzy_by_pos(bnd_lines, match_chrom):
    """For SV callers that do not report a MATEID such as sniffles2

    Here, the matching only works if the ALT allele reports
    the partner BND position exactly
    """

    bnd_lines = sorted(bnd_lines)
    # index all calls by their genomic location
    # and store list position for access
    breakend_locations = []
    for ln, _, bnd_line in bnd_lines:
        alt_spec = bnd_line[4]
        partner_pos = match_chrom.search(alt_spec)
        if partner_pos is None:
            call_id = BREAKENDS(bnd_line[0], bnd_line[1], None, None, ln)
        else:
            chrom2, pos2 = partner_pos.group().split(":")
            call_id = BREAKENDS(bnd_line[0], bnd_line[1], chrom2, pos2, ln)
        breakend_locations.append(call_id)

    # check for weird duplicates
    assert len(set(breakend_locations)) == len(breakend_locations)

    # this is a compromise / debug / workaround for sniffles,
    # which may have a bug to not properly report the positions
    # of the BND partners in a 1-to-1 manner. Likely can be dropped
    # if Sniffles has seen a fix
    # TODO - revise
    group_matches = set()
    # keep track of how many matches individual
    # BND calls attract to avoid comparing all vs all later
    match_counter = col.Counter()

    for bnd1, bnd2 in itt.combinations(breakend_locations, 2):
        assert bnd1.ln != bnd2.ln
        match12 = bnd1.chrom1 == bnd2.chrom2 and bnd1.pos1 == bnd2.pos2
        match21 = bnd2.chrom1 == bnd1.chrom2 and bnd2.pos1 == bnd1.pos2
        if match12 and match21:
            group_matches.add((bnd1.ln, bnd2.ln))
            match_counter[bnd1.ln] += 1
            match_counter[bnd2.ln] += 1

    if max(match_counter.values()) == 1:
        # each BND paired with at most one other

        bnd_by_line = dict(
            (ln, (rn, bnd_line)) for (ln, rn, bnd_line) in bnd_lines
        )

        collect_partners = col.defaultdict(set)
        group_members = col.defaultdict(set)
        # collect_partners[group_id].add((ln, rn, call_id, mate_id, bnd_line[0], tuple(bnd_line)))
        # group_members[group_id].add(call_id)

        group_counter = 1
        # NB: bnd_lines has been sorted above
        bnd_handled = set()
        for ln1, ln2 in sorted(group_matches):
            rn1, bnd1 = bnd_by_line[ln1]
            rn2, bnd2 = bnd_by_line[ln2]
            collect_partners[group_counter].add(
                (ln1, rn1, bnd1[2], bnd2[2], bnd1[0], tuple(bnd1))
            )
            group_members[group_counter].add(bnd1[2])
            collect_partners[group_counter].add(
                (ln2, rn2, bnd2[2], bnd1[2], bnd2[0], tuple(bnd2))
            )
            group_members[group_counter].add(bnd2[2])
            bnd_handled.add(ln1)
            bnd_handled.add(ln2)
            group_counter += 1


            #(ln, rn, call_id, mate_id, bnd_line[0], tuple(bnd_line)))
        for (ln, rn, bnd_line) in bnd_lines:
            if ln in bnd_handled:
                continue
            collect_partners[group_counter].add((ln, rn, bnd_line[2], None, bnd_line[0], tuple(bnd_line)))
            group_members[group_counter].add(bnd_line[2])
            group_counter += 1
    else:
        raise RuntimeError("Multi-paired breakends")

    return collect_partners, group_members


def derive_new_group_names(group_partners, group_members):


    bnd_index_info = []
    new_bnd_lines = []

    new_name_hashes = set()

    for group_id, partner_infos in group_partners.items():
        assert len(group_members[group_id]) == len(partner_infos)
        if len(partner_infos) == 1:
            # singleton
            bnd_type = "SNG"
            call_info = partner_infos.pop()
            assert call_info[3] is None
            bnd_vcf_line = "\t".join(call_info[-1])
            new_name = hl.md5((bnd_vcf_line).encode("utf-8")).hexdigest()
            assert new_name not in new_name_hashes
            new_name_hashes.add(new_name)
            new_name += bnd_type
            old_name = call_info[2]

            bnd_index_info.append((call_info[0], call_info[1], new_name, old_name))
            new_bnd_line = bnd_vcf_line.replace(old_name, new_name) + "\n"
            new_bnd_lines.append((call_info[0], new_bnd_line))
            continue

        group_hash_id = hl.md5("".join(sorted(group_members[group_id])).encode("utf-8")).hexdigest()
        assert group_hash_id not in new_name_hashes
        new_name_hashes.add(group_hash_id)

        group_chroms = set([t[4] for t in partner_infos])
        if len(group_chroms) > 1:
            bnd_type = "TLC"
        else:
            bnd_type = "ITR"

        group_hash_id += bnd_type
        total = len(partner_infos)

        new_member_ids = dict(
            (t[2], f"{group_hash_id}{i}v{total}") for i, t in enumerate(sorted(partner_infos), start=1)
        )

        for member_record in sorted(partner_infos):
            ln = member_record[0]  # line num
            rn = member_record[1]  # record num
            # sorted here is by input line number
            old_call_id = member_record[2]
            new_call_id = new_member_ids[old_call_id]
            old_mate_id = member_record[3]
            new_mate_id = new_member_ids[old_mate_id]

            bnd_vcf_line = "\t".join(member_record[-1]) + "\n"
            new_bnd_line = bnd_vcf_line.replace(
                old_call_id, new_call_id
            ).replace(old_mate_id, new_mate_id)

            new_bnd_lines.append((ln, new_bnd_line))
            bnd_index_info.append(((ln, rn, new_call_id, old_call_id)))

    return new_bnd_lines, bnd_index_info



def rename_bnd_partners(bnd_lines, has_mateid_field, match_chrom):

    if has_mateid_field:
        partner_infos, group_members = rename_exactly_by_mateid(bnd_lines)
    else:
        partner_infos, group_members = rename_fuzzy_by_pos(bnd_lines, match_chrom)

    new_bnd_lines, bnd_index_info = derive_new_group_names(partner_infos, group_members)

    if DEBUG:
        debug_msg(bnd_index_info)

    return new_bnd_lines, bnd_index_info


def main():

    args = parse_command_line()

    name_index = []
    # 2024-12-12
    # Some SV callers (for example, Sniffles2) do not report MATEIDS
    # for BND calls, which makes subsequent "joint" analysis ... difficult.
    # For Sniffles in particular, this is a known issue
    # (https://github.com/fritzsedlazeck/Sniffles/issues/121)
    # and, additionally, the type of the translocation may be
    # wrongly specified in ALT
    # https://github.com/fritzsedlazeck/Sniffles/issues/510

    # workaround / fix in general
    # If the input VCF is detected as being an SV callset
    # (look for SVTYPE field), we will buffer all lines
    # - hopefully just a few 10k of text ... - and determine
    # the respective BND pairings (NB: can be more than 2 partners!)
    # before renaming. This implies, though, that we also need to
    # buffer the header lines.
    sv_buffer = []
    bnd_line_buffer = []
    header_buffer = []

    is_sv_callset = None
    has_mateid_field = None

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
                header_buffer.append(line.strip())
                write_output.write(line)
                continue
            # at this point, the header buffer needs to be
            # non-empty, otherwise the VCF has not header
            if is_sv_callset is None:
                if not header_buffer:
                    raise RuntimeError("VCF input has no header")
                is_sv_callset, has_mateid_field, match_chrom = parse_vcf_header(header_buffer)
            record_num += 1

            if is_sv_callset:
                if "SVTYPE=BND" in line:
                    bnd_line_buffer.append((ln, record_num, line.strip().split()))
                else:
                    new_name = hl.md5(line.strip().encode("utf-8")).hexdigest()
                    columns = line.strip().split()
                    old_name = columns[2]
                    columns[2] = new_name
                    sv_buffer.append((ln, "\t".join(columns) + "\n"))
                    name_index.append((ln, record_num, new_name, old_name))
            else:
                new_name = hl.md5(line.strip().encode("utf-8")).hexdigest()
                columns = line.strip().split()
                old_name = columns[2]
                columns[2] = new_name
                write_output.write("\t".join(columns) + "\n")
                name_index.append((ln, record_num, new_name, old_name))

        if sv_buffer:
            assert is_sv_callset
            new_bnd_lines, bnd_index_entries = rename_bnd_partners(
                bnd_line_buffer, has_mateid_field, match_chrom
            )
            name_index.extend(bnd_index_entries)
            name_index = sorted(name_index)
            sv_buffer.extend(new_bnd_lines)
            # buffer now sorted by line numbers of input
            sv_buffer = sorted(sv_buffer)

            assert len(sv_buffer) == record_num

            [
                write_output.write(t[1]) for t in sv_buffer
            ]



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
