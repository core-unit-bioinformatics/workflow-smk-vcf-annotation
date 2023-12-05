"""
Use this module to list all includes
required for your pipeline - do not
add your pipeline-specific modules
to "commons/00_commons.smk"
"""

include: "00-prepare/settings.smk"
include: "00-prepare/sample_table.smk"
include: "00-prepare/ref_chroms.smk"

include: "10-norm/norm_index.smk"

include: "20-basic-filter/vemb_filters.smk"
include: "20-basic-filter/00_loc_pass.smk"
include: "20-basic-filter/10_genotype.smk"

include: "30-tabulate/pyutils.smk"
include: "30-tabulate/00_convert_vcf.smk"
include: "30-tabulate/10_vcf_tables.smk"

include: "40-merge/10_by_refpos.smk"
include: "40-merge/30_flatten_tables.smk"
include: "40-merge/50_group_calls.smk"
