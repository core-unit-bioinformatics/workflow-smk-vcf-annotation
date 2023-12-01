
rule split_vcf_tables:
    input:
        table = rules.tabulate_variant_size_distribution.output.table
    output:
        subset = DIR_PROC.joinpath(
            "30-tabulate", "10_vcf_tables",
            "00_split_chrom_group",
            "{sample}.{callset}.{ref}",
            "{sample}.{callset}.{ref}.basic-pass.{chrom}.{variant_group}.bed.gz"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "30-tabulate", "10_vcf_tables",
            "00_split_chrom_group",
            "{sample}.{callset}.{ref}",
            "{sample}.{callset}.{ref}.basic-pass.{chrom}.{variant_group}.split.rsrc"
        )
    shell:
        "zgrep -v chrom {input.table}"  # get rid of header
            " | "
        "{{ egrep \"^{wildcards.chrom}\\b\" || true; }}"  # do not fail if no calls on chrom
            " | "
        "{{ egrep \"\\s{wildcards.variant_group}\\s\" || true; }}"  # do not fail if no calls of type
            " | "
        "sort -V -k1 -k2n,3n"
            " | "
        "gzip -c > {output.subset}"


rule concat_all_vcf_subsets:
    input:
        subsets = get_callsets_by_ref
    output:
        concat = DIR_PROC.joinpath(
            "30-tabulate", "10_vcf_tables", "10_concat_by_ref",
            "{ref}", "{ref}.{chrom}.{variant_group}.bed.gz"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "30-tabulate", "10_vcf_tables", "10_concat_by_ref",
            "{ref}", "{ref}.{chrom}.{variant_group}.concat.rsrc"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "zcat {input.subsets} | sort -V -k1 -k2n,3n | gzip > {output.concat}"


# TODO add options to just concat by sample (= merge by caller)
