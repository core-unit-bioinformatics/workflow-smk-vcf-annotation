
rule create_grouped_long_call_table:
    input:
        tsv = expand(
            rules.flatten_merge_tables_long.output.multiples,
            chrom=config["reference_chromosomes"],
            allow_missing=True
        )
    output:
        tsv = DIR_RES.joinpath(
            "call_tables", "merged_groups", "{ref}",
            "{ref}.{variant_group}.by-group.tsv.gz"
        ),
    conda:
        DIR_ENVS.joinpath("vcftools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("merge_by_group")
    shell:
        "{params.script} --input-table {input.tsv} --output-table {output.tsv}"


rule count_sample_callset_subsets:
    input:
        tsv = rules.create_grouped_long_call_table.output.tsv,
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
        script=find_script("merge_by_group")
    shell:
        "{params.script} --grouped-calls {input.tsv} "
            "--out-subset-call-map {output.subset_map} "
            "--out-subsets {output.subset_counts}"


rule run_all_create_grouped_long_call_tables:
    input:
        tsv = expand(
            rules.create_grouped_long_call_table.output.tsv,
            ref=REFERENCE_GENOMES,
            variant_group=["SV"]
        ),
        counts = expand(
            rules.count_sample_callset_subsets.output.subset_counts,
            ref=REFERENCE_GENOMES,
            variant_group=["SV"]
        ),
        mapping = expand(
            rules.count_sample_callset_subsets.output.subset_map,
            ref=REFERENCE_GENOMES,
            variant_group=["SV"]
        ),
