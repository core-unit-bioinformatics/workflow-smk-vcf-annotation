
rule create_merged_group_indicator_table:
    input:
        tsv = expand(
            rules.flatten_merge_tables_long.output.multi_tsv,
            chrom=config["reference_chromosomes"],
            allow_missing=True
        )
    output:
        tsv = DIR_RES.joinpath(
            "callsets", "merged_groups", "{ref}", "tables",
            "{ref}.{variant_group}.group-indicator-table.tsv.gz"
        ),
    conda:
        DIR_ENVS.joinpath("vcftools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt * attempt,
        time_hrs=lambda wildcards, attempt: attempt * attempt
    params:
        script=find_script("merge_by_group")
    shell:
        "{params.script} --input-table {input.tsv} --output-table {output.tsv}"


rule dump_indicator_table_to_bedlike:
    input:
        tsv = rules.create_merged_group_indicator_table.output.tsv
    output:
        bed_like = DIR_RES.joinpath(
            "callsets", "merged_groups", "{ref}", "bed",
            "{ref}.{variant_group}.merged-groups.sample-sets.bed.gz"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt * attempt
    run:
        import pandas as pd
        select_columns = [
            "chrom", "start", "end", "group_id",
            "size", "vartype", "distinct_samples",
            "sample_set", "alt_allele_freq"
        ]
        df = pd.read_csv(input.tsv, sep="\t", header=0, usecols=select_columns)
        df.rename({"chrom": "#chrom"}, axis=1, inplace=True)
        df.to_csv(output.bed_like, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule count_sample_callset_subsets:
    input:
        tsv = rules.create_merged_group_indicator_table.output.tsv,
    output:
        subset_counts = DIR_RES.joinpath(
            "call_tables", "merged_groups", "{ref}",
            "{ref}.{variant_group}.subset-counts.tsv.gz"
        ),
        subset_map = DIR_RES.joinpath(
            "call_tables", "merged_groups", "{ref}",
            "{ref}.{variant_group}.subset-call-map.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("vcftools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        script=find_script("count_sample_sets")
    shell:
        "{params.script} --grouped-calls {input.tsv} "
            "--out-subset-call-map {output.subset_map} "
            "--out-subsets {output.subset_counts}"


rule run_all_create_merged_group_indicator_tables:
    input:
        tsv = expand(
            rules.create_merged_group_indicator_table.output.tsv,
            ref=REFERENCE_GENOMES,
            variant_group=["SV", "INDEL", "SNV"]
        ),
        bed = expand(
            rules.dump_indicator_table_to_bedlike.output.bed_like,
            ref=REFERENCE_GENOMES,
            variant_group=["SV", "INDEL", "SNV"]
        )
        # counts = expand(
        #     rules.count_sample_callset_subsets.output.subset_counts,
        #     ref=REFERENCE_GENOMES,
        #     variant_group=["SV"]
        # ),
        # mapping = expand(
        #     rules.count_sample_callset_subsets.output.subset_map,
        #     ref=REFERENCE_GENOMES,
        #     variant_group=["SV"]
        # ),
