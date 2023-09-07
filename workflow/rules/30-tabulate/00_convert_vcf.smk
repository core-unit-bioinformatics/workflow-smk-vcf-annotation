
rule tabulate_variant_size_distribution:
    input:
        vcf = rules.keep_variant_genotypes.output.vcf,
        tbi = rules.keep_variant_genotypes.output.tbi,
    output:
        table = DIR_PROC.joinpath(
            "30-tabulate", "00_convert_vcf",
            "{sample}.{callset}.{ref}.basic-pass.tsv.gz"
        ),
        stats = DIR_RES.joinpath(
            "statistics", "size_dist",
            "{sample}.{callset}.{ref}.basic-pass.stats.tsv"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "30-tabulate", "00_convert_vcf",
            "{sample}.{callset}.{ref}.basic-pass.tabconv.rsrc"
        ),
    wildcard_constraints:
        sample = CONSTRAINT_SAMPLES,
        callset = CONSTRAINT_CALLSETS,
        ref = CONSTRAINT_REFS
    conda:
        DIR_ENVS.joinpath("vcftools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt,
    params:
        script=find_script("tab_size_dist")
    shell:
        "{params.script} --vcf {input.vcf} --out-table {output.table} "
            "--out-size-dist {output.stats} --sample {wildcards.sample} "
            "--callset {wildcards.callset}"
