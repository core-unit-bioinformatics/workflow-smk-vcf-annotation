import pathlib

SAMPLES = set()
CALLSETS = set()
SAMPLE_CALLSET_MAP = dict()
SAMPLE_CALLSET_WILDCARDS = []
SAMPLE_CALLSET_REF_WILDCARDS = dict()

REFERENCE_GENOMES = dict()
for ref_id, ref_fasta in config["reference_genomes"].items():
    fasta_suffx = pathlib.Path(ref_fasta).suffix
    ref_infos = {
        "tag": ref_id,
        "fasta": ref_fasta,
        "faidx": pathlib.Path(ref_fasta).with_suffix(f"{fasta_suffx}.fai")
    }
    REFERENCE_GENOMES[ref_id] = ref_infos

input_path = pathlib.Path("/home/ebertp/work/projects/seq3105/data/callsets")

callsets = input_path.glob("*.vcf.gz")

for callset_vcf in callsets:
    assert callset_vcf.with_suffix(".gz.tbi").is_file()
    typed_sample = callset_vcf.name.split("_")[0]
    assert typed_sample[-1] in ["t", "b"]
    biotype = "tumor" if typed_sample[-1] == "t" else "blood"
    sample = typed_sample[:-1]
    callset_id = callset_vcf.name.rsplit(".", 2)[0]
    callset_id_parts = callset_id.split(".")
    callset_id = ".".join([callset_id_parts[1], callset_id_parts[3]])
    callset_ref = callset_id_parts[2]
    callset_id = f"{biotype}.{callset_id}"

    SAMPLES.add(sample)
    CALLSETS.add(callset_id)
    SAMPLE_CALLSET_MAP[(sample, callset_id, callset_ref)] = callset_vcf

    wc_comb = {
        "sample": sample,
        "callset": callset_id,
        "ref": callset_ref
    }
    SAMPLE_CALLSET_WILDCARDS.append(wc_comb)

SAMPLES = sorted(SAMPLES)
CALLSETS = sorted(CALLSETS)
REFERENCES = sorted(REFERENCE_GENOMES.keys())

CONSTRAINT_SAMPLES = "(" + "|".join(SAMPLES) + ")"
CONSTRAINT_CALLSETS = "(" + "|".join(CALLSETS) + ")"
CONSTRAINT_REFS = "(" + "|".join(REFERENCES) + ")"
CONSTRAINT_REFERENCES = CONSTRAINT_REFS

def get_sample_callset_wildcards(samples, callsets, references):
    return SAMPLE_CALLSET_WILDCARDS
