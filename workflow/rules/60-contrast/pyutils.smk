import collections
import itertools


def build_valid_contrast_group_combinations(*wildcards):

    flat_wildcards = collections.defaultdict(list)
    for wc_name, wc_value in itertools.chain.from_iterable(wildcards):
        flat_wildcards[wc_name].append(wc_value)

    all_refs = sorted(set(flat_wildcards["ref"]))
    all_variant_groups = sorted(set(flat_wildcards["variant_group"]))
    all_contrasts = sorted(set(flat_wildcards["contrast"]))

    valid_combinations = []
    for contrast in all_contrasts:
        group1_name = CONTRAST[contrast]["group1"]
        group2_name = CONTRAST[contrast]["group2"]

        for ref in all_refs:
            for vg in all_variant_groups:
                comb1 = {
                    "ref": ref,
                    "variant_group": vg,
                    "contrast": contrast,
                    "group_id": group1_name
                }
                valid_combinations.append(comb1)
                comb2 = dict(comb1)
                comb2["group_id"] = group2_name
                valid_combinations.append(comb2)
                comb3 = dict(comb2)
                comb3["group_id"] = "other"

    return sorted(valid_combinations)


def get_contrast_group_call_ids(wildcards):

    call_id_file = None
    contrast_spec = CONTRAST[wildcards.contrast]
    if wildcards.group_id == contrast_spec["group1"]:
        call_id_file = rules.get_grouped_multi_call_listing.output.calls_group1
    elif wildcards.group_id == contrast_spec["group2"]:
        call_id_file = rules.get_grouped_multi_call_listing.output.calls_group2
    elif wildcards.group_id == "other":
        call_id_file = rules.get_grouped_multi_call_listing.output.not_selected
    else:
        raise ValueError(f"Unknown contrast group: {wildcards}")
    assert call_id_file is not None
    return call_id_file
