import pathlib

import pandas

SAMPLES = set()
CALLSETS = set()
SAMPLE_CALLSET_MAP = dict()
SAMPLE_CALLSET_WILDCARDS = []


def parse_sample_sheet():

    sample_sheet = pandas.read_csv(
        SAMPLE_SHEET_PATH,
        sep="\t", header=0,
        comment="#"
    )

    global SAMPLES
    global CALLSETS
    global SAMPLE_CALLSET_MAP
    global SAMPLE_CALLSET_WILDCARDS

    for row in sample_sheet.itertuples():
        vcf_path = pathlib.Path(row.input_path).resolve(strict=True)
        tbi_path = pathlib.Path(row.input_path).with_suffix(".gz.tbi").resolve(strict=True)

        sample = row.sample
        callset_id = row.callset_id
        callset_ref = row.callset_ref

        SAMPLES.add(sample)
        CALLSETS.add(callset_id)
        SAMPLE_CALLSET_MAP[(sample, callset_id, callset_ref)] = vcf_path
        SAMPLE_CALLSET_MAP[(sample, callset_id, callset_ref, "vcf")] = vcf_path
        SAMPLE_CALLSET_MAP[(sample, callset_id, callset_ref, "idx")] = tbi_path

        wildcard_combination = {
            "sample": sample,
            "callset": callset_id,
            "ref": callset_ref
        }
        SAMPLE_CALLSET_WILDCARDS.append(wildcard_combination)

    return None


def build_constraint(wildcard_values):
    return "(" + "|".join(sorted(wildcard_values)) + ")"


def get_sample_callset_wildcards(samples, callsets, references):
    return SAMPLE_CALLSET_WILDCARDS


_ = parse_sample_sheet()


REFERENCE_MISMATCH = 0
for sc_key in SAMPLE_CALLSET_MAP.keys():
    # NB: key can have length 4 for vcf/idx info
    sample, callset_id, ref_id = sc_key[:3]
    if ref_id in REFERENCES:
        continue
    logerr(f"Annotated reference is not configured: {sample} / {callset_id} / {ref_id}")
    REFERENCE_MISMATCH += 1

if REFERENCE_MISMATCH > 0:
    logerr("Reference/callset mismatch - aborting.")
    logerr(f"{REFERENCE_GENOMES}")
    raise RuntimeError(
        f"Encountered {REFERENCE_MISMATCH} reference/callset mismatch error(s). "
        "See above for details."
    )


SAMPLES = sorted(SAMPLES)
CALLSETS = sorted(CALLSETS)

CONSTRAINT_SAMPLES = build_constraint(SAMPLES)
CONSTRAINT_CALLSETS = build_constraint(CALLSETS)
