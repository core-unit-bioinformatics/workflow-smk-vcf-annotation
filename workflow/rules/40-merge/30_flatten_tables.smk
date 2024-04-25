
rule flatten_merge_tables_short:
    input:
        tsv = rules.merge_identical_short_by_refpos.output.merged
    output:
        single_tsv = DIR_RES.joinpath(
            "callsets", "singletons", "{ref}", "tables",
            "{ref}.{chrom}.{variant_group}.singletons.flat.tsv.gz"
        ),
        multi_tsv = DIR_RES.joinpath(
            "callsets", "multicalls", "{ref}", "tables",
            "{ref}.{chrom}.{variant_group}.multicalls.flat.tsv.gz"
        ),
        count_stats = DIR_RES.joinpath(
            "statistics", "counts", "{ref}", "by_chrom",
            "{ref}.{chrom}.{variant_group}.flat-tables.count-stats.tsv"
        ),
    wildcard_constraints:
        variant_group="(INDEL|SNV)"
    conda:
        DIR_ENVS.joinpath("vcftools.yaml")
    params:
        script=find_script("flatten_strict")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt,
        time_hrs=lambda wildcards, attempt: attempt
    shell:
        "{params.script} --merged-table {input.tsv} --variant-group {wildcards.variant_group} "
            "--singletons {output.single_tsv} --multicalls {output.multi_tsv} "
            "--count-stats {output.count_stats}"


rule flatten_merge_tables_long:
    input:
        infos = rules.concat_vcf_subsets_by_ref_chrom.output.concat,
        mrg = rules.merge_proximal_long_by_refpos.output.merged,
    output:
        single_tsv = DIR_RES.joinpath(
            "callsets", "singletons", "{ref}", "tables",
            "{ref}.{chrom}.{variant_group}.singletons.flat.tsv.gz"
        ),
        multi_tsv = DIR_RES.joinpath(
            "callsets", "multicalls", "{ref}", "tables",
            "{ref}.{chrom}.{variant_group}.multicalls.flat.tsv.gz"
        ),
        count_stats = DIR_RES.joinpath(
            "statistics", "counts", "{ref}", "by_chrom",
            "{ref}.{chrom}.{variant_group}.flat-tables.count-stats.tsv"
        ),
        desc_stats = DIR_PROC.joinpath(
            "statistics", "descriptive", "{ref}", "by_chrom",
            "{ref}.{chrom}.{variant_group}.flat-tables.desc-stats.tsv"
        ),
    wildcard_constraints:
        variant_group="SV"
    conda:
        DIR_ENVS.joinpath("vcftools.yaml")
    params:
        script=find_script("sv_merging")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt,
        time_hrs=lambda wildcards, attempt: attempt
    shell:
        "{params.script} --sv-table {input.infos} --sv-merge {input.mrg} "
            "--singletons {output.single_tsv} --multicalls {output.multi_tsv} "
            "--count-stats {output.count_stats} --desc-stats {output.desc_stats} "


rule run_all_flatten_merge_tables_short:
    input:
        singles = expand(
            rules.flatten_merge_tables_short.output.single_tsv,
            ref=REFERENCES,
            chrom=config["reference_chromosomes"],
            variant_group=["INDEL", "SNV"]

        ),
        multis = expand(
            rules.flatten_merge_tables_short.output.multi_tsv,
            ref=REFERENCES,
            chrom=config["reference_chromosomes"],
            variant_group=["INDEL", "SNV"]
        ),
        stats = expand(
            rules.flatten_merge_tables_short.output.count_stats,
            ref=REFERENCES,
            chrom=config["reference_chromosomes"],
            variant_group=["INDEL", "SNV"]
        )


rule run_all_flatten_merge_tables_long:
    input:
        singles = expand(
            rules.flatten_merge_tables_long.output.single_tsv,
            ref=REFERENCES,
            chrom=config["reference_chromosomes"],
            variant_group=["SV"]
        ),
        multis = expand(
            rules.flatten_merge_tables_long.output.multi_tsv,
            ref=REFERENCES,
            chrom=config["reference_chromosomes"],
            variant_group=["SV"]
        ),
        count_stats = expand(
            rules.flatten_merge_tables_long.output.count_stats,
            ref=REFERENCES,
            chrom=config["reference_chromosomes"],
            variant_group=["SV"]
        ),
        desc_stats = expand(
            rules.flatten_merge_tables_long.output.desc_stats,
            ref=REFERENCES,
            chrom=config["reference_chromosomes"],
            variant_group=["SV"]
        )

