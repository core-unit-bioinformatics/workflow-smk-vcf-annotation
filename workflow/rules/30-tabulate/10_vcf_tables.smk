
rule split_vcf_tables_by_chrom_group:
    input:
        table = rules.convert_vcf_to_table.output.table
    output:
        subset = DIR_PROC.joinpath(
            "30-tabulate", "10_vcf_tables",
            "00_split_chrom_group",
            "{sample}.{callset}.{ref}",
            "{sample}.{callset}.{ref}.basic-pass.{chrom}.{variant_group}.tsv.gz"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "30-tabulate", "10_vcf_tables",
            "00_split_chrom_group",
            "{sample}.{callset}.{ref}",
            "{sample}.{callset}.{ref}.basic-pass.{chrom}.{variant_group}.split.rsrc"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    shell:
        "zgrep -e \"^chrom\\s\" {input.table} | gzip > {output.subset}"
            " && "
        "zgrep -v chrom {input.table}"  # skip over header
            " | "
        "{{ egrep \"^{wildcards.chrom}\\b\" || true; }}"  # do not fail if no calls on chrom
            " | "
        "{{ egrep \"\\s{wildcards.variant_group}\\s\" || true; }}"  # do not fail if no calls of type
            " | "
        "sort -V -k1 -k2n,3n"
            " | "
        "gzip >> {output.subset}"


rule concat_vcf_subsets_by_ref_chrom:
    input:
        subsets = get_callsets_by_ref
    output:
        concat = DIR_PROC.joinpath(
            "30-tabulate", "10_vcf_tables", "10_concat_by_ref",
            "{ref}", "{ref}.{chrom}.{variant_group}.tsv.gz"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "30-tabulate", "10_vcf_tables", "10_concat_by_ref",
            "{ref}", "{ref}.{chrom}.{variant_group}.concat.rsrc"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "zgrep -e \"^chrom\\s\" {input.subsets[0]} | gzip > {output.concat}"
            " && "
        "zcat {input.subsets}"
            " | "
        "egrep -v \"^chrom\""  # skip over header lines
            " | "
        "sort -V -k1 -k2n,3n | gzip >> {output.concat}"


rule concat_chrom_callsets_by_ref:
    """This rule is just a convenience exit point
    of the workflow that essentially "copies" all
    concatenated sample-/callset-level calls into
    the result folder for (potential) later inspection.
    """
    input:
        tables = expand(
            rules.concat_vcf_subsets_by_ref_chrom.output.concat,
            chrom=config["reference_chromosomes"],
            allow_missing=True
        )
    output:
        table = DIR_RES.joinpath(
            "callsets", "concat_by_ref", "tables",
            "{ref}", "{ref}.{variant_group}.concat-calls.tsv.gz"
        )
    resources:
        mem_mb=lambda wildcards, attempt: {
            "SV": 2048,
            "SNV": 8192,
            "INDEL": 4096
        }[wildcards.variant_group] * attempt
    shell:
        "zgrep -e \"^chrom\\s\" {input.tables[0]} | gzip > {output.table}"
            " && "
        "zcat {input.tables}"
            " | "
        "egrep -v \"^chrom\""  # skip over header lines
            " | "
        "sort -V -k1 -k2n,3n | gzip >> {output.table}"


rule run_all_concat_vcf_subsets_by_ref:
    """TODO: reference chromosomes should become a workflow variable
    """
    input:
        tables = expand(
            rules.concat_vcf_subsets_by_ref_chrom.output.concat,
            ref=REFERENCES,
            chrom=config["reference_chromosomes"],
            variant_group=["SNV", "INDEL", "SV"]
        ),
        res_table = expand(
            rules.concat_chrom_callsets_by_ref.output.table,
            ref=REFERENCES,
            variant_group=["SNV", "INDEL", "SV"]
        )


# TODO add options to just concat by sample (= merge by caller)
