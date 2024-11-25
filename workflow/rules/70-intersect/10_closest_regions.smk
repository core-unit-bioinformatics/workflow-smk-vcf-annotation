
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
            ANNOTATIONS[wildcards.ref][wildcards.annotation]
        ),
        chrom_list = rules.write_ref_chrom_lists.output.regexp
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
    resources:
        mem_mb=lambda wildcards, attempt: (8192 if wildcards.variant_group == "SNV" else 2048) * attempt
    run:
        import gzip
        import pathlib as pl
        import pandas as pd

        callset_header = gzip.open(input.call_bed, "rt").readline().strip().split()
        assert callset_header[0].startswith("#")

        ref_suffix = pl.Path(input.ref_bed).suffix
        if ref_suffix == ".bed":
            open_ref, open_mode = open, "r"
        elif ref_suffix == ".gz":
            open_ref, open_mode = gzip.open, "rt"
        else:
            logerr(f"Unexpected file format: {ref_suffix}")
            raise ValueError(f"Cannot read file format: {ref_suffix}")

        ann_header = open_ref(input.ref_bed, open_mode).readline().strip().split()
        if ann_header[0].startswith("#"):
            ann_header[0] = f"chrom"
        try:
            _ = int(ann_header[1])
            _ = int(ann_header[2])
            # columns 1 and 2 are integers --- not a header line
            if len(ann_header) == 4:
                ann_header = ["chrom", "start", "end", "name"]
            else:
                logerr(f"Annotation file has no valid header: {input.ref_bed}")
                raise RuntimeError(f"No valid header line in {input.ref_bed}")
        except ValueError:
            # columns 1 and 2 are not integers --- likely valid header line
            pass
        ann_header = [f"{column}_{wildcards.annotation}" for column in ann_header]

        header_intersect = set(callset_header).intersection(set(ann_header))
        if len(header_intersect) > 0:
            logerr(f"Key / column label collision: {header_intersect}")
            raise RuntimeError(f"Incompatible headers: {input.call_bed} / {input.ref_bed}")

        input_table_header = callset_header + ann_header + [f"distance_{wildcards.annotation}"]
        df = pd.read_csv(input.tsv, sep="\t", header=None, names=input_table_header)

        df.to_csv(output.bed_like, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule merge_and_filter_annotations:
    input:
        positives = lambda wildcards: expand(
            rules.reheader_intersect_tables.output.bed_like,
            annotation=ANNOTATION_FILTERS[wildcards.ref]["positive"],
            allow_missing=True
        ),
        negatives = lambda wildcards: expand(
            rules.reheader_intersect_tables.output.bed_like,
            annotation=ANNOTATION_FILTERS[wildcards.ref]["negative"],
            allow_missing=True
        )
    output:
        bed_like = DIR_RES.joinpath(
            "annotations", "merged_filtered", "closest_region",
            "{ref}.{variant_group}.{callset_type}.closest.{annotation}.bed.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt
    params:
        script=find_script("join_ann_tables.py")
    shell:
        "{params.script} --positive {input.positives} --negative {input.negative} --output {output.bed_like}"


rule run_all_find_closest_annotated_region:
    # TODO
    # this must be turned into a configurable pull
    input:
        tables_hg38_sv = expand(
            rules.reheader_intersect_tables.output.bed_like,
            ref=["hg38"],
            annotation=sorted(ANNOTATIONS["hg38"].keys()),
            variant_group=["SV"],
            callset_type=[
                "groupcalls", "bycatch-25-75", "bycatch-10-90"
            ] + CONTRAST_CALLSET_LABELS
        ),
        tables_hg38_short = expand(
            rules.reheader_intersect_tables.output.bed_like,
            ref=["hg38"],
            annotation=sorted(ANNOTATIONS["hg38"].keys()),
            variant_group=["SNV", "INDEL"],
            callset_type=[
                "groupcalls"
            ] + CONTRAST_CALLSET_LABELS
        ),
        tables_t2t_sv = expand(
            rules.reheader_intersect_tables.output.bed_like,
            ref=["t2tv2"],
            annotation=sorted(ANNOTATIONS["t2tv2"].keys()),
            variant_group=["SV"],
            callset_type=[
                "groupcalls", "bycatch-25-75", "bycatch-10-90"
            ] + CONTRAST_CALLSET_LABELS
        ),
        tables_t2t_short = expand(
            rules.reheader_intersect_tables.output.bed_like,
            ref=["t2tv2"],
            annotation=sorted(ANNOTATIONS["t2tv2"].keys()),
            variant_group=["SNV", "INDEL"],
            callset_type=[
                "groupcalls"
            ] + CONTRAST_CALLSET_LABELS
        ),
        merged_filtered = expand(
            rules.merge_and_filter_annotations.output.bed_like,
            ref=["hg38"],
            variant_group=["SV"],
            callset_type=["groupcalls"] + CONTRAST_CALLSET_LABELS
        )
