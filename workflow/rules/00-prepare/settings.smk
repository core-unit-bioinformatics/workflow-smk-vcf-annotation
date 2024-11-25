
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

KNOWN_VARIANT_GROUPS = ["SV", "SNV", "INDEL"]
CONSTRAINT_VAR_GROUPS = "(" + "|".join(KNOWN_VARIANT_GROUPS) + ")"


# collect annotations from config
_ANNOTATIONS = config.get("annotations", None)
if _ANNOTATIONS is None:
    ANNOTATIONS = dict()
else:
    ANNOTATIONS = dict()
    for ref_genome, annotations in _ANNOTATIONS.items():
        assert ref_genome in REFERENCES
        if ref_genome not in ANNOTATIONS:
            ANNOTATIONS[ref_genome] = dict()
        for label, ann_file in annotations.items():
            assert label not in ANNOTATIONS[ref_genome]
            ANNOTATIONS[ref_genome][label] = ann_file


# collect annotation filters from config
_ANNOTATION_FILTERS = config.get("annotation_filters", None)
if _ANNOTATION_FILTERS is None:
    ANNOTATION_FILTERS = dict()
else:
    ANNOTATION_FILTERS = dict()
    for ref_genome, filter_lists in ANNOTATION_FILTERS.items():
        pos_list = filter_lists.get("positive", [])
        neg_list = filter_lists.get("negative", [])
        ANNOTATION_FILTERS[ref_genome] = {
            "positive": pos_list,
            "negative": neg_list
        }
