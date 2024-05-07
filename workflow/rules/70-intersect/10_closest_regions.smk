
rule norm_enforce_sort_order_callset:
    input:
        callset = lambda wildcards: get_bedlike_callset(wildcards.callset_type)
    output:
        callset = temp(
            DIR_PROC.joinpath(
                "temp", "enforce_sort_order",
                "{ref}.{variant_group}.{callset_type}.bed.gz"
            )
        )
    params:
        grep = lambda wildcards, input: select_grep_cmd(input.callset)
    shell:
        "{params.grep} -e -v \"^#\" {input.callset}"
            " | "
        "sort -V -k1,1 -k2,2n -k3,3n"
            " | "
        "gzip > {output.callset}"


rule norm_enforce_sort_order_annotation:
    input:
        bedlike = lambda wildcards: DIR_GLOBAL_REF.joinpath(
            config["annotations"][wildcards.ref][wildcards.annotation]
        )
    output:
        bedlike = temp(
            DIR_PROC.joinpath(
                "temp", "enforce_sort_order",
                "{ref}.{annotation}.bed.gz"
        ))
    params:
        grep = lambda wildcards, input: select_grep_cmd(input.bedlike)
    shell:
        "{params.grep} -e -v \"^#\" {input.bedlike}"
            " | "
        "sort -V -k1,1 -k2,2n -k3,3n"
            " | "
        "gzip > {output.bedlike}"


rule find_closest_annotated_region:
    """Since bedtools is very stubborn and as flexible as
    an anvil when it comes to the chromosome sort order in
    the input files, this rule relies on presorted temporary
    input that is created as a preprocessing step.
    """
    input:
        calls = rules.norm_enforce_sort_order_callset.output.callset,
        ann = rules.norm_enforce_sort_order_annotation.output.bedlike,
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
