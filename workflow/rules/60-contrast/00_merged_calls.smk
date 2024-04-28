

rule get_grouped_multi_call_listing:
    """For grouped multi calls, it is fairly straightforward
    to disjoin the calls by group given the sample_set column
    subset mapping.
    Why the generic identifiers group1 and group2?
    """
    input:
        merged_calls = rules.create_merged_group_indicator_table.output.tsv
    output:
        calls_group1 = DIR_RES.joinpath(
            "contrast", "{ref}",
            "{ref}.{variant_group}.contrast-{contrast}-group1.group-call-ids.txt"
        ),
        calls_group2 = DIR_RES.joinpath(
            "contrast", "{ref}",
            "{ref}.{variant_group}.contrast-{contrast}-group2.group-call-ids.txt"
        ),
        not_selected = DIR_RES.joinpath(
            "contrast", "{ref}",
            "{ref}.{variant_group}.contrast-{contrast}-other.group-call-ids.txt"
        ),
        contrast_info = DIR_RES.joinpath(
            "contrast", "{ref}",
            "{ref}.{variant_group}.contrast-{contrast}.group-info.txt"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    run:
        import pandas as pd
        df = pd.read_csv(input.merged_calls, sep="\t", header=0, usecols=["group_id", "sample_set"])

        def check_sample_group_match(subset_samples, group1_samples, group2_samples):
            """Criterion check:
            - at least 1 of the positive samples is in the subset
            - none of the negative samples is in the subset
            ---> count subset as hit/match for positive samples
            """
            samples = subset_samples.split(",")
            group1_match = any(group1_sample in samples for group1_sample in group1_samples)
            group2_match = any(group2_sample in samples for group2_sample in group2_samples)
            if group1_match and group2_match:
                # matching samples in both groups
                group_label = 0
            elif group1_match:  # and not group2
                group_label = 1
            elif group2_match:  # and not group2
                group_label = 2
            else:
                # no sample matches for either group
                group_label = 0
            return group_label

        group1_samples = CONTRAST[wildcards.contrast]["samples1"]
        group2_samples = CONTRAST[wildcards.contrast]["samples2"]

        group_labels = df["subset"].apply(check_sample_group_match, args=(group1_samples, group2_samples))
        df["group_label"] = group_labels

        select_group1 = df["group_label"] == 1
        assert select_group1.any()
        group1_calls = set().union(
            *df.loc[select_group1, "call_group_ids"].apply(lambda call_ids: set(call_ids.split(",")))
        )
        select_group2 = df["group_label"] == 2
        assert select_group2.any()
        group2_calls = set().union(
            *df.loc[select_group2, "call_group_ids"].apply(lambda call_ids: set(call_ids.split(",")))
        )

        select_other = df["group_label"] == 0
        assert select_other.any()
        other_calls = set().union(
                *df.loc[select_other, "call_group_ids"].apply(lambda call_ids: set(call_ids.split(",")))
        )

        # fully disjoin call sets - see docstring for explanation
        g1g2_intersect = group1_calls.intersection(group2_calls)
        g1oth_intersect = group1_calls.intersection(other_calls)
        g2oth_intersect = group2_calls.intersection(other_calls)
        group1_calls = group1_calls - g1g2_intersect - g1oth_intersect
        group2_calls = group2_calls - g1g2_intersect - g2oth_intersect
        other_calls = other_calls.union(g1_g2_intersect)

        assert len(group1_calls.intersection(group2_calls)) == 0
        assert len(group1_calls.intersection(other_calls)) == 0
        assert len(group2_calls.intersection(other_calls)) == 0
        assert len(group1_calls) + len(group2_calls) + len(other_calls) == df.shape[0]

        callsets = [group1_calls, group2_calls, other_calls]
        outfiles = [output.calls_group1, calls_group2, not_selected]

        for calls, outfile in zip(callsets, outfiles):
            with open(outfile, "w") as dump:
                _ = dump.write("\n".join(sorted(calls)) + "\n")

        with open(output.contrast_info, "w") as dump:
            _ = dump.write(f"contrast_name\t{wildcards.contrast}\n")
            _ = dump.write(f"contrast_group1\t{CONTRAST[wildcards.contrast]['group1']}\n")
            _ = dump.write(f"contrast_samples1\t{','.join(CONTRAST[wildcards.contrast]['samples1'])}\n")
            _ = dump.write(f"contrast_group2\t{CONTRAST[wildcards.contrast]['group2']}\n")
            _ = dump.write(f"contrast_samples2\t{','.join(CONTRAST[wildcards.contrast]['samples2'])}\n")
    # END OF RUN BLOCK


rule extract_contrast_merged_group_calls:
    input:
        table = rules.create_merged_group_indicator_table.output.tsv,
        call_ids = get_contrast_group_call_ids
    output:
        table = DIR_RES.joinpath(
            "callsets", "contrasts", "{ref}", "tables",
            "{ref}.{variant_group}.contrast-{contrast}-{group_id}.group-indicator-table.tsv.gz"
        )
    wildcard_constraints:
        group_id="(group1|group2|other)"
    shell:
        "zcat {input.table} | egrep \"^chrom\" | gzip > {output.table}"
            " && "
        "zgrep -F -f {input.call_ids} {input.table} | gzip >> {output.table}"


rule dump_contrast_indicator_table_to_bedlike:
    input:
        tsv = rules.extract_contrast_merged_group_calls.output.table
    output:
        bed_like = DIR_RES.joinpath(
            "callsets", "contrasts", "{ref}", "bed",
            "{ref}.{variant_group}.contrast-{contrast}-{group_id}.merged-groups.sample-sets.bed.gz"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt * attempt
    run:
        import pandas as pd
        select_columns = [
            "chrom", "start", "end", "group_id",
            "size", "vartype", "distinct_samples", "sample_set"
        ]
        df = pd.read_csv(input.tsv, sep="\t", header=0, usecols=select_columns)
        df.rename({"chrom": "#chrom"}, axis=1, inplace=True)
        df.to_csv(output.bed_like, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_build_contrast_id_lists:
    input:
        bed_like = expand(
            rules.dump_contrast_indicator_table_to_bedlike.output.bed_like,
            build_valid_contrast_group_combinations,
            ref=REFERENCE_GENOMES,
            variant_group=["SV", "INDEL"],
            contrast=list(CONTRAST.keys()),
            group_id=CONTRAST_GROUPS1 + CONTRAST_GROUPS2 + ["other"]
        )
