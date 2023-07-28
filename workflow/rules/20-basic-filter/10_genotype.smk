
rule keep_variant_genotypes:
    input:
        vcf = DIR_PROC.joinpath(
            "10-filter/{sample}/00_loc_pass/{sample}.{callset}.loc-pass.vcf.gz"
        ),
        tbi = DIR_PROC.joinpath(
            "10-filter/{sample}/00_loc_pass/{sample}.{callset}.loc-pass.vcf.gz.tbi"
        ),
    output:
        vcf = DIR_PROC.joinpath(
            "10-filter/{sample}/10_genotype/{sample}.{callset}.gt-ok.vcf.gz"
        ),
        tbi = DIR_PROC.joinpath(
            "10-filter/{sample}/10_genotype/{sample}.{callset}.gt-ok.vcf.gz.tbi"
        ),
    log:
        DIR_LOG.joinpath(
            "10-filter/10_genotype/{sample}.{callset}.keep-gt.vembrane.log"
        )
    conda:
        DIR_ENVS.joinpath("vcftools.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    params:
        filter_expr = get_genotype_filter()
    shell:
        "vembrane filter --output-fmt vcf"
            " {params.filter_expr} {input.vcf} 2> {log}"
            " | "
        "bcftools sort --max-mem 2G"
            " | "
        "bgzip --stdout > {output.vcf}"
            " && "
        "tabix -p vcf {output.vcf}"


rule exclude_ref_missing_genotypes:
    input:
        vcf = DIR_PROC.joinpath(
            "10-filter/{sample}/00_loc_pass/{sample}.{callset}.loc-pass.vcf.gz"
        ),
        tbi = DIR_PROC.joinpath(
            "10-filter/{sample}/00_loc_pass/{sample}.{callset}.loc-pass.vcf.gz.tbi"
        ),
    output:
        vcf = DIR_PROC.joinpath(
            "10-filter/{sample}/10_genotype/excluded",
            "{sample}.{callset}.excl-gt.vcf.gz"
        ),
        tbi = DIR_PROC.joinpath(
            "10-filter/{sample}/10_genotype/excluded",
            "{sample}.{callset}.excl-gt.vcf.gz.tbi"
        ),
    log:
        DIR_LOG.joinpath(
            "10-filter/10_genotype/{sample}.{callset}.excl-gt.vembrane.log"
        )
    conda:
        DIR_ENVS.joinpath("vcftools.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    params:
        filter_expr = get_genotype_filter(exclude=True)
    shell:
        "vembrane filter --output-fmt vcf"
            " {params.filter_expr} {input.vcf} 2> {log}"
            " | "
        "bcftools sort --max-mem 2G"
            " | "
        "bgzip --stdout > {output.vcf}"
            " && "
        "tabix -p vcf {output.vcf}"


rule run_basic_filter_genotype:
    input:
        vcf_keep = expand(
            rules.keep_variant_genotypes.output.vcf,
            get_sample_callset_wildcards,
            sample=SAMPLES,
            callset=CALLSETS
        ),
        vcf_excl = expand(
            rules.exclude_ref_missing_genotypes.output.vcf,
            get_sample_callset_wildcards,
            sample=SAMPLES,
            callset=CALLSETS
        ),
