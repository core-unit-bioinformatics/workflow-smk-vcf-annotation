
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
