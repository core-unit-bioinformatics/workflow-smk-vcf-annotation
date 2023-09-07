
rule merge_snv_by_ref:
    """Since SNVs do not need
    to be checked for overlap criteria,
    we can just merge identical variants
    (by location in this step).
    """
    input:
        concat = rules.concat_all_vcf_subsets.output.concat
    output:
        merged = DIR_PROC.joinpath(
            "40-merge", "snv_tables", "merge_by_ref",
            "{ref}.{chrom}.{variant_group}.merged.bed.gz"
        )
    wildcard_constraints:
        variant_group="SNV"
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "bedtools merge -delim \"|\" -c 6,7,8,20 -o collapse -i {input.concat}"
            " | "
        "pigz -p {threads} > {output.merged}"


rule run_merge_snv_by_ref:
    input:
        merged = expand(
            rules.merge_snv_by_ref.output.merged,
            ref=["t2tv2", "hg38"],
            chrom=config["reference_chromosomes"],
            variant_group=["SNV"]
        )
