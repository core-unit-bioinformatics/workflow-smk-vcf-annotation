
rule merge_identical_short_by_refpos:
    """For InDels and SNVs, apply a conservative strategy
    (by default) and merge only identical calls.
    This is driven by the position in the reference
    (bedtools merge command) and then refined during
    a post-processing step: indels of different lengths
    would be split in different groups or become
    singletons.
    """
    input:
        concat = rules.concat_all_vcf_subsets.output.concat
    output:
        merged = DIR_PROC.joinpath(
            "40-merge", "10_by_refpos", "{ref}",
            "{ref}.{chrom}.{variant_group}.merged.tsv.gz"
        )
    wildcard_constraints:
        variant_group="(SNV|INDEL)"
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        select_columns=lambda wildcards: {
            "SNV": "6,7,8,16,20",
            "INDEL": "4,5,6,7,8,16,20"
        }[wildcards.variant_group]
    shell:
        "bedtools merge -delim \"|\" -d -1 -c {params.select_columns} "
        "-o collapse -i {input.concat}"
            " | "
        "pigz -p {threads} > {output.merged}"


rule run_merge_identical_short_by_refpos:
    input:
        merged = expand(
            rules.merge_identical_short_by_refpos.output.merged,
            ref=REFERENCES,
            chrom=config["reference_chromosomes"],
            variant_group=["INDEL", "SNV"]
        )
