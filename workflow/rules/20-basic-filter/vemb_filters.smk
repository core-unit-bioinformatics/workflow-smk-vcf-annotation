

def get_loc_pass_filter(exclude=False):
    """
    Vembrane filter expression that filters
    for reference chromosomes (location)
    and all-pass FILTER conditions (or ".")
    """

    keep_filter = (
        "'"
        "CHROM in AUX[\"chroms\"] and "
        "all(f in [\"PASS\", \".\"] for f in FILTER)"
        "'"
    )

    exclude_filter = (
        "'"
        "CHROM not in AUX[\"chroms\"] or "
        "all(f not in [\"PASS\", \".\"] for f in FILTER)"
        "'"
    )

    select_filter = exclude_filter if exclude else keep_filter

    return select_filter


def get_genotype_filter(exclude=False):
    """
    Vembrane filter expression that filters
    for uninteresting genotypes, i.e.,
    all calls that do not involve a variant type
    or where ALT contains more than one unresolved
    sequence symbol (N)
    - the latter criterion applies to cuteSV callsets,
    but must not check for SUM(N in ALT) == 0 because
    sniffles sets ALT alleles to N
    """

    keep_filter = (
        "'"
        "count_any_var() > 0 "
        "and "
        "sum(c == \"N\" for c in ALT.upper()) < 2 "
        "and "
        "sum(c == \"N\" for c in REF.upper()) < 2"
        "'"
    )

    exclude_filter = (
        "'"
        "count_any_var() < 1 "
        "or "
        "sum(c == \"N\" for c in ALT.upper()) > 1 "
        "or "
        "sum(c == \"N\" for c in REF.upper()) > 1"
        "'"
    )

    select_filter = exclude_filter if exclude else keep_filter

    return select_filter
