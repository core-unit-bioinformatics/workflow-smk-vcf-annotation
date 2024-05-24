
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
    wildcard_constraints:
        ref=CONSTRAINT_REFS,
        variant_group=CONSTRAINT_VAR_GROUPS
    params:
        grep = lambda wildcards, input: select_grep_cmd(input.callset)
    shell:
        "{params.grep} -E -v \"^#\" {input.callset}"
            " | "
        "sort -V -k1,1 -k2,2n -k3,3n"
            " | "
        "gzip > {output.callset}"


rule norm_enforce_sort_order_annotation:
    input:
        bedlike = lambda wildcards: DIR_GLOBAL_REF.joinpath(
            config["annotations"][wildcards.ref][wildcards.annotation]
        ),
        chrom_list = rules.write_ref_chrom_lists.output.listing
    output:
        bedlike = temp(
            DIR_PROC.joinpath(
                "temp", "enforce_sort_order",
                "{ref}.{annotation}.bed.gz"
        ))
    params:
        grep = lambda wildcards, input: select_grep_cmd(input.bedlike)
    shell:
        "{params.grep} -E -v \"^#\" {input.bedlike}"
            " | "
        "grep -w -f {input.chrom_list}"
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


rule reheader_intersect_tables:
    input:
        tsv = rules.find_closest_annotated_region.output.tsv,
        ref_bed = lambda wildcards: DIR_GLOBAL_REF.joinpath(
            config["annotations"][wildcards.ref][wildcards.annotation]
        ),
        call_bed = lambda wildcards: get_bedlike_callset(wildcards.callset_type)
    output:
        bed_like = DIR_RES.joinpath(
            "annotations", "separate", "closest_region",
            "{ref}.{variant_group}.{callset_type}.closest.{annotation}.bed.gz"
        )
    wildcard_constraints:
        ref=CONSTRAINT_REFS,
        variant_group=CONSTRAINT_VAR_GROUPS
    run:
        import gzip
        import pathlib as pl
        import pandas as pd

        callset_header = gzip.open(input.call_bed, "rt").readline().strip().split()
        assert callset_header[0].startswith("#")
        last_callset_column = len(callset_header)

        ref_suffix = pl.Path(input.ref_bed).suffix
        if ref_suffix == ".bed":
            open_ref, open_mode = open, "r"
        elif ref_suffix == ".gz":
            open_ref, open_mode = gzip.open, "rt"
        else:
            logerr(f"Unexpected file format: {ref_suffix}")
            raise ValueError(f"Cannot read file format: {ref_suffix}")

        ref_header = open_ref(input.ref_bed, open_mode).readline().strip().split()
        if not ref_header[0].startswith("#"):
            if len(ref_header) == 4:
                ref_header = ["chrom2", "start2", "end2", "ref_name"]
            else:
                logerr(f"Reference file w/o header: {input.ref_bed}")
                raise ValueError(f"Malformed reference file (>4 columns, no header): {input.ref_bed}")
        else:
            # this would commonly be something like #chrom
            ref_header[0] = "chrom2"
            ref_header = list(map(str.lower, ref_header))

        header_intersect = set(callset_header).intersection(set(ref_header))
        if len(header_intersect) > 0:
            logerr(f"Callset and reference header not disjoint: {header_intersect}")
            raise ValueError(
                "Callset and reference header must be disjoined:\n"
                f" {input.call_bed} / {input.ref_bed} / {header_intersect}"
            )

        input_table_header = callset_header + ref_header + ["distance"]
        df = pd.read_csv(input.tsv, sep="\t", header=None, names=input_table_header)
        assert "annotation" not in input_table_header
        df.insert(last_callset_column, "annotation", wildcards.annotation)

        df.to_csv(output.bed_like, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_find_closest_annotated_region:
    # TODO
    # this must be turned into a configurable pull
    input:
        tables = expand(
            rules.reheader_intersect_tables.output.bed_like,
            ref=["hg38"],
            annotation=["genes", "hgsvc2", "bands", "ogm", "cosmic"],
            variant_group=["SV"],
            callset_type=[
                "groupcalls", "bycatch-25-75", "bycatch-10-90"
            ] + CONTRAST_CALLSET_LABELS
        ),
        tables2 = expand(
            rules.reheader_intersect_tables.output.bed_like,
            ref=["t2tv2"],
            annotation=["bands", "newseq", "genes", "uniq"],
            variant_group=["SV"],
            callset_type=[
                "groupcalls", "bycatch-25-75", "bycatch-10-90"
            ] + CONTRAST_CALLSET_LABELS
        )
