
def get_bedlike_callset(callset_type):

    if callset_type == "groupcalls":
        # 40-merge::50_group_calls.smk
        callset_bed = rules.dump_indicator_table_to_bedlike.output.bed_like
    elif callset_type == "multicalls":
        raise NotImplementedError()
    elif callset_type == "singletons":
        raise NotImplementedError()
    elif "bycatch" in callset_type:
        # 50-augment::20_save_singletons.smk
        # NB: "bycatch" callsets are only defined for SV calls at the moment
        bracket = callset_type.split("-", 1)[-1]
        callset_bed = expand(
            rules.concat_long_singletons_bycatch.output.bed,
            bracket=bracket,
            allow_missing=True
        )
    else:
        raise ValueError(f"Unknown callset type: {callset_type}")

    return callset_bed


def get_contrast_bedlike_callset():
    raise NotImplementedError()
