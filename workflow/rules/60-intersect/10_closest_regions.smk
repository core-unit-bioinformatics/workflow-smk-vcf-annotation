
rule reduce_table_to_bed:
    input:
        tsv = rules.create_grouped_long_call_table.output.tsv
    output:
        bed = DIR_RES.joinpath(
            "call_tables", "merged_groups", "{ref}",
            "{ref}.{variant_group}.by-group.bed.gz"
        )
    shell:
        "zcat {input.tsv} | tail -n +2 | cut -f 1-6 | gzip > {output.bed}"


rule find_closest_annotated_region:
    input:
        calls = rules.reduce_table_to_bed.output.bed,
        ann = lambda wildcards: DIR_GLOBAL_REF.joinpath(
            config["annotations"][wildcards.ref][wildcards.annotation]
        )
    output:
        tsv = DIR_PROC.joinpath(
            "60-intersect", "closest_regions",
            "{ref}.{variant_group}.grouped-calls.{annotation}.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "bedtools closest -d -a {input.calls} -b {input.ann}"
            " | "
        "gzip > {output.tsv}"


rule run_all_find_closest_annotated_region:
    input:
        tables = expand(
            rules.find_closest_annotated_region.output.tsv,
            ref=["hg38"],
            annotation=["genes", "hgsvc2"],
            variant_group=["SV"]
        )
