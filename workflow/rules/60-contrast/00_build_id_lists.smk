

rule get_grouped_multi_call_listing:
    """For grouped multi calls, it is fairly straightforward
    to disjoin the calls by group given the sample/callset
    subset mapping. NB:
    (1) the subset map first needs to be reduced to just sample
    subsets (and not callsets) to properly define the 'not_selected' output.
    (2) the subset map lists shared calls for all possible subsets, i.e.
    the sample subset (A,B) has a list of shared calls, but such a call could
    also be shared with another sample, e.g. subset (A,C). Hence, to properly
    disjoin the calls between the two groups, subtracting the set intersections
    at the end of the code block is essential.
    """
    input:
        call_map = rules.count_sample_callset_subsets.output.subset_map
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
        )
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    run:
        import pandas as pd
        df = pd.read_csv(input.call_map, sep="\t", header=0)
        select_sample_subsets = df["subset"].apply(lambda subset: any(s in subset for s in SAMPLES))
        df = df.loc[select_sample_subsets, :].copy()

        def check_sample_group_match(subset_samples, group_samples):
            samples = subset_samples.split(",")
            group_match = all(s in group_samples for s in samples)
            return group_match

        group1_samples = CONTRAST[wildcards.contrast]["samples1"]
        select_group1 = df["subset"].apply(check_sample_group_match, args=(group1_samples,))

        group2_samples = CONTRAST[wildcards.contrast]["samples2"]
        select_group2 = df["subset"].apply(check_sample_group_match, args=(group2_samples,))
        select_other = ~(select_group1 or select_group2)

        group1_calls = set().union(
            *df.loc[select_group1, "call_group_ids"].apply(lambda call_ids: set(call_ids.split(",")))
        )
        group2_calls = set().union(
            *df.loc[select_group2, "call_group_ids"].apply(lambda call_ids: set(call_ids.split(",")))
        )
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

        callsets = [group1_calls, group2_calls, other_calls]
        outfiles = [output.calls_group1, calls_group2, not_selected]

        for calls, outfile in zip(callsets, outfiles):
            with open(outfile, "w") as dump:
                _ = dump.write("\n".join(sorted(calls)) + "\n")
    # END OF RUN BLOCK


rule run_all_build_contrast_id_lists:
    input:
        listings = expand(
            rules.get_grouped_multi_call_listing.output,
            ref=REFERENCE_GENOMES,
            variant_group=["SV"],
            contrast=list(CONTRAST.keys())
        )
