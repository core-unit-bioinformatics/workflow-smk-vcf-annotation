
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
        "zcat {input.table} | head -1 | gzip > {output.subset}"
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
        "zcat {input.subsets[0]} | head -1 | gzip > {output.concat}"
            " && "
        "zgrep -v chrom {input.subsets} | sort -V -k1 -k2n,3n | gzip >> {output.concat}"


rule run_all_concat_vcf_subsets_by_ref:
    """TODO: reference chromosomes should become a workflow variable
    """
    input:
        tables = expand(
            rules.concat_vcf_subsets_by_ref_chrom.output.concat,
            ref=REFERENCES,
            chrom=config["reference_chromosomes"],
            variant_group=["SNV", "INDEL", "SV"]
        )


# TODO add options to just concat by sample (= merge by caller)
