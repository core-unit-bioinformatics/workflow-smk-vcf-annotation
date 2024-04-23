
VCF_NORM_REF_ACTION = config.get("vcf_norm_ref_action", "s")
assert VCF_NORM_REF_ACTION in ["s", "e", "x", "w"]

REFERENCE_GENOMES = dict()
for ref_id, ref_fasta in config["reference_genomes"].items():
    fasta_suffx = pathlib.Path(ref_fasta).suffix
    ref_infos = {
        "tag": ref_id,
        "fasta": ref_fasta,
        "faidx": pathlib.Path(ref_fasta).with_suffix(f"{fasta_suffx}.fai")
    }
    REFERENCE_GENOMES[ref_id] = ref_infos

REFERENCES = sorted(REFERENCE_GENOMES.keys())
CONSTRAINT_REFS = "(" + "|".join(REFERENCES) + ")"
CONSTRAINT_REFERENCES = CONSTRAINT_REFS

#########################################
# Process user-defined contrasts
# Contrasts specify pairwise comparisons
# (samples or groups of samples)
# to auto-create callsets of interest,
# i.e. calls in one group vs the other.
#########################################

CONTRAST = dict()
CONTRAST_GROUPS1 = []
CONTRAST_GROUPS2 = []
_CONTRAST_STRUCT = config.get("contrast", None)
if _CONTRAST_STRUCT is not None:
    for contrast_name, contrast_groups in _CONTRAST_STRUCT.items():
        if isinstance(contrast_groups, dict):
            if len(contrast_groups) != 2:
                err_msg = (
                    "Currently, only group- or sample-pairwise contrasts are supported. "
                    f"You did not specify exactly two groups: {list(contrast_groups.keys())}"
                )
                logerr(err_msg)
                raise ValueError(err_msg)

            contrast_spec = {
                "name": contrast_name
            }
            n = 1
            for group, samples in contrast_groups.items():
                if samples[0] == "auto":
                    if n != 2:
                        err_msg = (
                            "Automatic sample collection can only be set for second group. "
                            f"Error in contrast: {contrast_name}"
                        )
                        logerr(err_msg)
                        raise ValueError(err_msg)
                    samples1 = contrast_spec["samples1"]
                    samples2 = [s for s in SAMPLES if s not in samples1]
                    contrast_spec["group2"] = group
                    contrast_spec["samples2"] = samples2
                    continue
                contrast_spec[f"group{n}"] = group
                contrast_spec[f"samples{n}"] = samples
                n += 1
            assert n == 2
            CONTRAST[contrast_name] = contrast_spec

        elif isinstance(contrast_groups, list):
            if len(contrast_groups) != 2 or not all(s in SAMPLES for s in contrast_groups):
                err_msg = (
                    "Detected simple list as contrast group, which is only allowed "
                    "to specify a sample-pair contrast. The list must have length two "
                    "and all items in the list must be valid sample identifiers. "
                    f"Your list is: {contrast_groups}"
                )
                logerr(err_msg)
                raise ValueError(err_msg)
            contrast_spec = {
                "name": contrast_name,
                "group1": contrast_groups[0],
                "samples1": contrast_groups[0],
                "group2": contrast_groups[1],
                "samples2": contrast_groups[1]
            }
            assert contrast_name not in CONTRAST, f"Duplicated contrast name: {contrast_name}"
            CONTRAST[contrast_name] = contrast_spec
        else:
            raise ValueError(f"Cannot parse contrast group definition: {contrast_groups}")

    # populate lists CONTRAST_GROUPS1 and CONTRAST_GROUPS2
    # to be used as wildcard replacements
    for contrast, contrast_spec in CONTRAST.items():
        CONTRAST_GROUPS1.append(contrast_spec["group1"])
        CONTRAST_GROUPS2.append(contrast_spec["group2"])
