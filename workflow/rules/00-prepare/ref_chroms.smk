
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


localrules: write_genome_size_file
rule write_genome_size_file:
    input:
        idx = lambda wildcards: DIR_GLOBAL_REF.joinpath(REFERENCE_GENOMES[wildcards.ref]["faidx"])
    output:
        gsize = DIR_LOCAL_REF.joinpath("{ref}.size")
    run:
        chrom_sizes = []
        with open(input.idx, "r") as listing:
            for line in listing:
                chrom, size = line.split()[:2]
                chrom_sizes.append((chrom, int(size)))
        chrom_sizes = sorted(chrom_sizes, key=lambda t: t[1], reverse=True)
        with open(output.gsize, "w") as listing:
            for chrom, size in chrom_sizes:
                _ = listing.write(f"{chrom}\t{size}\n")
    # END OF RUN BLOCK
