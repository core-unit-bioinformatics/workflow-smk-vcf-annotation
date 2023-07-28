
rule normalize_variant_calls:
    """bcftools 1.17 usage note:
    the options "--rm-dup" and "--multiallelics"
    are mutually exclusive (undocumented in CLI help)
    """
    input:
        vcf = lambda wildcards: SAMPLE_CALLSET_MAP[(wildcards.sample, wildcards.callset)],
        tbi = lambda wildcards: pathlib.Path(
            SAMPLE_CALLSET_MAP[(wildcards.sample, wildcards.callset)]
        ).with_suffix(".gz.tbi"),
        ref = lambda wildcards: DIR_GLOBAL_REF.joinpath(
            REF_GENOME_FOR_CALLSET[(wildcards.sample, wildcards.callset)]["fasta"]
        ),
    output:
        vcf = DIR_PROC.joinpath(
            "10-norm/{sample}.{callset}.norm-idx.vcf.gz"
        ),
        tbi = DIR_PROC.joinpath(
            "10-norm/{sample}.{callset}.norm-idx.vcf.gz.tbi"
        ),
        tsv = DIR_RES.joinpath(
            "aux", "record_index", "{sample}.{callset}.idx.tsv.gz"
        ),
        json = DIR_RES.joinpath(
            "aux", "record_index", "{sample}.{callset}.idx.json"
        )
    log:
        DIR_LOG.joinpath(
            "10-norm", "{sample}.{callset}.norm-idx.bcftools.log"
        )
    wildcard_constraints:
        sample = CONSTRAINT_SAMPLES,
        callset = CONSTRAINT_CALLSETS
    conda:
        DIR_ENVS.joinpath("vcftools.yaml")
    threads: CPU_LOW
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    params:
        script = find_script("build_vcf_name_index"),
        action = VCF_NORM_REF_ACTION
    shell:
        "bcftools norm --threads {threads} --atomize --atom-overlaps . "
            "--check-ref {params.action} --fasta-ref {input.ref} "
            "--multiallelics -any --multi-overlaps . "
            "--output-type v {input.vcf} 2> {log}"
            " | "
        "bcftools norm --threads {threads} --rm-dup exact "
            "--output-type v /dev/stdin 2>> {log}"
            " | "
        "bcftools sort --max-mem 2G "
            " | "
        "{params.script} --name-idx {output.tsv} --idx-info {output.json} "
        "--tags sample:{wildcards.sample} callset:{wildcards.callset}"
            " | "
        "bgzip --stdout > {output.vcf}"
            " && "
        "tabix -p vcf {output.vcf}"


rule run_normalize_all_callsets:
    input:
        vcf = expand(
            rules.normalize_variant_calls.output.vcf,
            get_sample_callset_wildcards,
            sample=SAMPLES,
            callset=CALLSETS
        )
