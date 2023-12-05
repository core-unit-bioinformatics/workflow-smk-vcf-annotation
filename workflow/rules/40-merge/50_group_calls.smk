
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
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("merge_by_group")
    shell:
        "{params.script} --input-table {input.tsv} --output-table {output.tsv}"


rule run_all_create_grouped_long_call_tables:
    input:
        tsv = expand(
            rules.create_grouped_long_call_table.output.tsv,
            ref=REFERENCE_GENOMES,
            variant_group=["SV"]
        )
