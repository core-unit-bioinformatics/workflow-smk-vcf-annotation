
localrules: write_ref_chrom_lists
rule write_ref_chrom_lists:
    input:
        ref = lambda wildcards: DIR_GLOBAL_REF.joinpath(REFERENCE_GENOMES[wildcards.ref]["fasta"]),
        idx = lambda wildcards: DIR_GLOBAL_REF.joinpath(REFERENCE_GENOMES[wildcards.ref]["faidx"]),
    output:
        listing = DIR_PROC.joinpath(
            "00-prepare/ref_chroms/{ref}.keep-chroms.txt"
        ),
    wildcard_constraints:
        ref=CONSTRAINT_REFERENCES
    run:
        ref_chroms = config["reference_chromosomes"]
        # TODO support var sets
        with open(output.listing, "w") as dump:
            _ = dump.write("\n".join(ref_chroms) + "\n")
    # END OF RUN BLOCK



