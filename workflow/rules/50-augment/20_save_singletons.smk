
rule select_singletons_in_bracket_long:
    """NB: the 'read_header' input is just >some<
    table to get the full call table header w/o
    resorting to ugly hacks such as realized in the
    'sv_merging.py' script (TODO - change that as well)

    call_table = 30-tabulate::10_vcf_tables::RULE
    singletons = 40-merge::30_flatten_tables::RULE
    """
    input:
        read_header = expand(
            rules.tabulate_variant_size_distribution.output.table,
            sample=SAMPLE_CALLSET_WILDCARDS[0]["sample"],
            callset=SAMPLE_CALLSET_WILDCARDS[0]["callset"],
            ref=SAMPLE_CALLSET_WILDCARDS[0]["ref"]
        ),
        call_table = rules.concat_all_vcf_subsets.output.concat,
        singletons = rules.flatten_merge_tables_long.output.singletons,
        stats_bracket = rules.determine_plausibility_thresholds_long.output.single_select_bracket
    output:
        sng_select_tsv = DIR_RES.joinpath(
            "call_tables", "singletons", "{ref}",
            "{ref}.{chrom}.{variant_group}.singletons.{bracket}.tsv.gz"
        ),
        sng_select_bed = DIR_RES.joinpath(
            "call_tables", "singletons", "{ref}",
            "{ref}.{chrom}.{variant_group}.singletons.{bracket}.bed.gz"
        )
    wildcard_constraints:
        variant_group="SV"
    conda:
        DIR_ENVS.joinpath("pyutils.yaml")
    params:
        script=find_script("get_singletons")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    shell:
        "{params.script} --read-header {input.read_header} --call-table {input.call_table} "
        "--singletons {input.singletons} --stats-bracket {input.stats_bracket} "
        "--out-table {output.sng_select_tsv} --out-bed {output.sng_select_bed}"


rule run_all_save_singletons_long:
    input:
        tsv = expand(
            rules.select_singletons_in_bracket_long.output.sng_select_tsv,
            ref=REFERENCES,
            chrom=config["reference_chromosomes"],
            variant_group=["SV"],
            bracket=["25-75", "10-90"]
        )
