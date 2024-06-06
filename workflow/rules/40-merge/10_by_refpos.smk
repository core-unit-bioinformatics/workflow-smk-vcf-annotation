
rule prep_merge_table_to_bed:
    """
    """
    input:
        table = rules.concat_vcf_subsets_by_ref_chrom.output.concat
    output:
        bed_like = temp(
            DIR_PROC.joinpath(
                "40-merge", "10_by_refpos", "00_tmp_bed",
                "{ref}", "{ref}.{chrom}.{variant_group}.concat-calls.bed.gz"
            )
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd

        select_columns = {
            "SV": ["name"],
            "SNV": ["name", "sample", "callset", "ref_allele_repr", "alt_allele_repr", "alt_allele_freq"],
            "INDEL": ["vartype", "size", "name", "sample", "callset", "ref_allele_repr", "alt_allele_repr", "alt_allele_freq"]
        }
        read_columns = ["chrom", "start", "end"] + select_columns[wildcards.variant_group]
        df = pd.read_csv(input.table, sep="\t", header=0, usecols=read_columns)
        df.sort_values(["start", "end"], inplace=True)
        df.rename({"chrom": "#chrom"}, axis=1, inplace=True)
        df.to_csv(output.bed_like, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


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
        concat = rules.prep_merge_table_to_bed.output.bed_like
    output:
        merged = DIR_PROC.joinpath(
            "40-merge", "10_by_refpos", "{ref}",
            "{ref}.{chrom}.{variant_group}.pos-merged.tsv.gz"
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
            "SNV": "4,5,6,7,8,9",
            "INDEL": "4,5,6,7,8,9,10,11"
        }[wildcards.variant_group]  # column def: see rule prep_merge_table_to_bed
    shell:
        "bedtools merge -header -delim \"|\" -d -1 -c {params.select_columns} "
        "-o collapse -i {input.concat}"
            " | "
        "pigz -p {threads} > {output.merged}"


rule merge_proximal_long_by_refpos:
    """For SV calls, perfect identity of two calls
    is too idealistic as merge criterion. Hence,
    apply a distance-based cutoff to first group
    calls in a simple bedtools merge step, and then
    examine the groups to finally decide if calls
    should be merged.
    TODO --- properly parameterize distance cutoff
    and merge into rule above / create single rule
    """
    input:
        concat = rules.prep_merge_table_to_bed.output.bed_like
    output:
        merged = DIR_PROC.joinpath(
            "40-merge", "10_by_refpos", "{ref}",
            "{ref}.{chrom}.{variant_group}.pos-merged.tsv.gz"
        )
    wildcard_constraints:
        variant_group="SV"
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        select_columns="4",  # column def: see rule prep_merge_table_to_bed
        dist_cutoff=config.get("sv_proximal_dist_cutoff", 200)
    shell:
        "bedtools merge -header -delim \"|\" -d {params.dist_cutoff} -c {params.select_columns} "
        "-o collapse -i {input.concat}"
            " | "
        "pigz > {output.merged}"


rule run_all_merge_identical_short_by_refpos:
    input:
        merged = expand(
            rules.merge_identical_short_by_refpos.output.merged,
            ref=REFERENCES,
            chrom=config["reference_chromosomes"],
            variant_group=["INDEL", "SNV"]
        )


rule run_all_merge_proximal_long_by_refpos:
    input:
        merged = expand(
            rules.merge_proximal_long_by_refpos.output.merged,
            ref=REFERENCES,
            chrom=config["reference_chromosomes"],
            variant_group=["SV"]
        )
