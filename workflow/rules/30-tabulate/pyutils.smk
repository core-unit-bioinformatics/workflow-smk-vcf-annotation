
# TODO: candidate for inclusion in workflow template
class SafeFormattingDict(dict):
    def __missing__(self, key):
        return "{" + key + "}"


def get_callsets_by_ref(wildcards):

    input_template = rules.split_vcf_tables.output.subset
    select_ref = wildcards.ref

    callset_for_ref = []
    for wildcard_combination in SAMPLE_CALLSET_WILDCARDS:
        if wildcard_combination["ref"] != select_ref:
            continue
        fmt_input = input_template.format_map(
            SafeFormattingDict(
                sample=wildcard_combination["sample"],
                callset=wildcard_combination["callset"],
            )
        )
        callset_for_ref.append(fmt_input)
    return callset_for_ref
