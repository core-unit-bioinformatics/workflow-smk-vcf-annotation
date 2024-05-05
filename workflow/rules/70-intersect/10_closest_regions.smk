
rule find_closest_annotated_region:
    input:
        calls = lambda wildcards: get_bedlike_callset(wildcards.callset_type),
        ann = lambda wildcards: DIR_GLOBAL_REF.joinpath(
            config["annotations"][wildcards.ref][wildcards.annotation]
        )
    output:
        tsv = DIR_PROC.joinpath(
            "70-intersect", "closest_regions",
            "{ref}.{variant_group}.{callset_type}.closest.{annotation}.tsv.gz"
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
    # TODO
    # this must be turned into a configurable pull
    input:
        tables = expand(
            rules.find_closest_annotated_region.output.tsv,
            ref=["hg38"],
            annotation=["genes", "hgsvc2", "bands", "ogm", "cosmic"],
            variant_group=["SV"],
            callset_type=["groupcalls", "bycatch-25-75", "bycatch-10-90"]
        ),
        tables2 = expand(
            rules.find_closest_annotated_region.output.tsv,
            ref=["t2tv2"],
            annotation=["bands", "newseq", "genes", "uniq"],
            variant_group=["SV"],
            callset_type=["groupcalls", "bycatch-25-75", "bycatch-10-90"]
        )
