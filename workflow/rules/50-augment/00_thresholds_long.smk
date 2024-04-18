
localrules: determine_plausibility_thresholds_long
rule determine_plausibility_thresholds_long:
    """Although highest confidence should be attributed
    to calls that are confirmed in other samples or
    at least by other callers, in highly variable samples
    (e.g., cancer), singletons should not just be discarded.
    This rule determines 'plausibility' thresholds to
    select singleton calls that fall within a conservative
    range for various statistics that are computed on the
    basis of the merged 'multi' calls
    """
    input:
        desc_stats = rules.flatten_merge_tables_long.output.desc_stats
    output:
        single_select_bracket = DIR_RES.joinpath(
            "statistics", "single_select_bracket",
            "{ref}", "{ref}.{chrom}.{variant_group}.sng-select.{bracket}.tsv"
        )
    params:
        t_low = lambda wildcards: int(wildcards.bracket.split("-")[0]),
        t_high = lambda wildcards: int(wildcards.bracket.split("-")[1])
    run:
        import pandas as pd
        _this_fun = "50-augment::00_thresholds_long::determine_plausibility_thresholds_long"

        df = pd.read_csv(input.desc_stats, sep="\t", header=0)
        bracket_low = f"pctile_{params.t_low}"
        bracket_high = f"pctile_{params.t_high}"
        check_low_exists = bracket_low in df["statistic"].values
        check_high_exists = bracket_high in df["statistic"].values
        if not (check_low_exists and check_high_exists):
            err_msg = (
                f"{_this_fun}\n"
                "Cannot determine plausibility thresholds; at least one "
                "percentile bracket value is not part of this statistics "
                f"file: {input.desc_stats}\n"
                f"Low percentile value: {params.t_low} - exists? {check_low_exists}\n"
                f"High percentile value: {params.t_high} - exists? {check_high_exists}\n"
            )
            logerr(err_msg)
            raise ValueError(err_msg)

        select_occurrence = df["occurrence"] == "multiple"
        select_lower = df["statistic"] == bracket_low
        select_upper = df["statistic"] == bracket_high
        selector = select_occurrence & select_lower & select_upper
        sub = df.loc[selector, :].copy()
        sub.to_csv(output.single_select_bracket, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


