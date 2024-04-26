
rule select_singletons_in_bracket_long:
    """NB: the call table must be the complete/pre-merged
    table because singletons that were already dropped
    due to the initial clustering by distance will not be
    be part of the post-merged singletons table.

    All singletons:
        call_table = 30-tabulate::10_vcf_tables::RULE

    Singletons that were not merged w/ other calls
    but that are close to a merged group call:
        singletons = 40-merge::30_flatten_tables::RULE
    """
    input:
        call_table = rules.concat_vcf_subsets_by_ref_chrom.output.concat,
        singletons = rules.flatten_merge_tables_long.output.single_tsv,
        stats_bracket = rules.determine_plausibility_thresholds_long.output.single_select_bracket
    output:
        sng_select_tsv = DIR_PROC.joinpath(
            "50-augment", "20_save_singletons", "{ref}", "by_chrom",
            "{ref}.{chrom}.{variant_group}.singletons.bycatch-{bracket}.tsv.gz"
        ),
        sng_select_bed = DIR_PROC.joinpath(
            "50-augment", "20_save_singletons", "{ref}", "by_chrom",
            "{ref}.{chrom}.{variant_group}.singletons.bycatch-{bracket}.bed.gz"
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
        "{params.script} --call-table {input.call_table} "
        "--singletons {input.singletons} --stats-bracket {input.stats_bracket} "
        "--out-table {output.sng_select_tsv} --out-bed {output.sng_select_bed}"


rule concat_long_singletons_bycatch:
    input:
        tables = expand(
            rules.select_singletons_in_bracket_long.output.sng_select_tsv,
            chrom=config["reference_chromosomes"],
            allow_missing=True
        ),
        bed = expand(
            rules.select_singletons_in_bracket_long.output.sng_select_bed,
            chrom=config["reference_chromosomes"],
            allow_missing=True
        ),
    output:
        table = DIR_RES.joinpath(
            "callsets", "singletons", "{ref}", "tables",
            "{ref}.{variant_group}.singletons.bycatch-{bracket}.tsv.gz"
        ),
        bed = DIR_RES.joinpath(
            "callsets", "singletons", "{ref}", "bed",
            "{ref}.{variant_group}.singletons.bycatch-{bracket}.tsv.gz"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt * attempt
    run:
        import pandas as pd
        all_inputs = [sorted(input.tables), sorted(input.bed)]
        all_outputs = [output.table, output.bed]

        for infiles, outfile in zip(all_inputs, all_outputs):
            concat = []
            for infile in infiles:
                df = pd.read_csv(infile, sep="\t", header=0)
                df.sort_values(["start", "end"], inplace=True)
                concat.append(df)
            concat = pd.concat(concat, axis=0, ignore_index=False)
            concat.to_csv(outfile, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_save_singletons_long:
    input:
        tsv = expand(
            rules.concat_long_singletons_bycatch.output,
            ref=REFERENCES,
            variant_group=["SV"],
            bracket=["25-75", "10-90"]
        )
