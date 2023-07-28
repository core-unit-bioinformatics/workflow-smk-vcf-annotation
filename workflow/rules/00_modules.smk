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
