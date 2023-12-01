
rule keep_all_pass_in_regions:
    input:
        vcf = DIR_PROC.joinpath(
            "10-norm/{sample}.{callset}.{ref}.norm-idx.vcf.gz"
        ),
        tbi = DIR_PROC.joinpath(
            "10-norm/{sample}.{callset}.{ref}.norm-idx.vcf.gz.tbi"
        ),
        keep_chroms = DIR_PROC.joinpath(
            "00-prepare/ref_chroms/{ref}.keep-chroms.txt"
        ),
    output:
        vcf = DIR_PROC.joinpath(
            "20-basic-filter/{sample}/00_loc_pass/{sample}.{callset}.{ref}.loc-pass.vcf.gz"
        ),
        tbi = DIR_PROC.joinpath(
            "20-basic-filter/{sample}/00_loc_pass/{sample}.{callset}.{ref}.loc-pass.vcf.gz.tbi"
        ),
    log:
        DIR_LOG.joinpath(
            "20-basic-filter/00_loc_pass/{sample}.{callset}.{ref}.keep-loc-pass.vembrane.log"
        )
    wildcard_constraints:
        sample = CONSTRAINT_SAMPLES,
        callset = CONSTRAINT_CALLSETS,
        ref = CONSTRAINT_REFS
    conda:
        DIR_ENVS.joinpath("vcftools.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    params:
        filter_expr = get_loc_pass_filter()
    shell:
        "vembrane filter --output-fmt vcf --aux chroms={input.keep_chroms}"
            " {params.filter_expr} {input.vcf} 2> {log}"
            " | "
        "bcftools sort --max-mem 2G"
            " | "
        "bgzip --stdout > {output.vcf}"
            " && "
        "tabix -p vcf {output.vcf}"


rule exclude_no_pass_outside_regions:
    input:
        vcf = DIR_PROC.joinpath(
            "10-norm/{sample}.{callset}.{ref}.norm-idx.vcf.gz"
        ),
        tbi = DIR_PROC.joinpath(
            "10-norm/{sample}.{callset}.{ref}.norm-idx.vcf.gz.tbi"
        ),
        keep_chroms = DIR_PROC.joinpath(
            "00-prepare/ref_chroms/{ref}.keep-chroms.txt"
        ),
    output:
        vcf = DIR_PROC.joinpath(
            "20-basic-filter/{sample}/00_loc_pass/excluded",
            "{sample}.{callset}.{ref}.excl-loc-pass.vcf.gz"
        ),
        tbi = DIR_PROC.joinpath(
            "20-basic-filter/{sample}/00_loc_pass/excluded",
            "{sample}.{callset}.{ref}.excl-loc-pass.vcf.gz.tbi"
        ),
    log:
        DIR_LOG.joinpath(
            "20-basic-filter/00_loc_pass/{sample}.{callset}.{ref}.excl-loc-pass.vembrane.log"
        )
    conda:
        DIR_ENVS.joinpath("vcftools.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    params:
        filter_expr = get_loc_pass_filter(exclude=True)
    shell:
        "vembrane filter --output-fmt vcf --aux chroms={input.keep_chroms}"
            " {params.filter_expr} {input.vcf} 2> {log}"
            " | "
        "bcftools sort --max-mem 2G"
            " | "
        "bgzip --stdout > {output.vcf}"
            " && "
        "tabix -p vcf {output.vcf}"


rule run_basic_filter_loc_pass:
    input:
        vcf_keep = expand(
            rules.keep_all_pass_in_regions.output.vcf,
            get_sample_callset_wildcards,
            sample=SAMPLES,
            callset=CALLSETS,
            ref=REFERENCES
        ),
        vcf_excl = expand(
            rules.exclude_no_pass_outside_regions.output.vcf,
            get_sample_callset_wildcards,
            sample=SAMPLES,
            callset=CALLSETS,
            ref=REFERENCES
        ),
