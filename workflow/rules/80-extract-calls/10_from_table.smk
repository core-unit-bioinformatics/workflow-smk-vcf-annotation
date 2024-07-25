

localrules: dump_group_call_ids
rule dump_group_call_ids:
    output:
        listing = DIR_PROC.joinpath(
            "80-extract-calls", "user_input",
            "{subset_name}.group-call-ids.txt"
        )
    run:
        import io
        import re
        # group and call IDs are MD5 checksums
        match_id = re.compile("[a-z0-9]{32}", flags=re.IGNORECASE)
        subset_lists = config.get("extract_group_calls", None)
        if subset_lists is None:
            err_msg = (
                "Rule 80-extract-calls::dump_group_call_ids triggered "
                "but no subset(s) of group calls specified in config with "
                "key 'extract_group_calls'. Aborting ..."
            )
            logerr(err_msg)
            raise ValueError(err_msg)

        this_subset = subset_lists.get(wildcards.subset_name, None)
        if this_subset is None:
            err_msg = (
                "Rule 80-extract-calls::dump_group_call_ids triggered "
                "but the user-specified subset of group calls with the "
                f"name '{wildcards.subset_name}' is not part of the "
                "respective config file. Aborting ..."
            )
            logerr(err_msg)
            raise ValueError(err_msg)

        out_buffer = io.StringIO()
        for call_id in this_subset:
            if match_id.match(call_id) is None:
                err_msg = (
                    "Rule 80-extract-calls::dump_group_call_ids encountered "
                    f"a malformed group call ID ({call_id}) as part of the "
                    f"subset listing {wildcards.subset_name}."
                )
                logerr(err_msg)
                raise ValueError(err_msg)
            out_buffer.write(f"{call_id}\n")

        with open(output.listing, "w") as dump:
            _ = dump.write(out_buffer.getvalue())
    # END OF RUN BLOCK


rule fetch_group_call_entries:
    input:
        tables = expand(
            output.create_grouped_long_call_table.output.tsv,
            ref=REFERENCE_GENOMES,
            variant_group=["SV", "INDEL", "SNV"]
        ),
        id_list = rules.dump_group_call_ids.output.listing
    output:
        subset = DIR_RES.joinpath(
            "user_subsets", "{subset_name}.tsv.gz"
        )
    shell:
        "zgrep -f {input.id_list} {input.tables}"
            " | "
        "gzip > {output.subset}"


rule run_all_extract_group_calls:
    input:
        subsets = expand(
            rule.fetch_group_call_entries.output.subset,
            subset_name=EXTRACT_GROUP_CALL_IDS
        )
