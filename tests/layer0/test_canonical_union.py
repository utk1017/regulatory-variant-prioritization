"""Locks in Problem #3 from the review: the original notebook restricted VEP parsing to a single
hardcoded canonical transcript, which caused a real 5' UTR variant to be mislabeled because the
canonical transcript's TSS definition was offset from a valid alternative transcript's. Layer 0
must union consequence terms across all of the target gene's transcripts, not just canonical."""
from regvar.layer0.vep import parse_vep_results
from tests.conftest import make_vep_result

TARGET_GENE = "ENSG00000244734"
CANONICAL_TX = "ENST00000335295"


def test_noncanonical_transcript_term_is_not_lost():
    """A term present only on a non-canonical (but same-gene) transcript must still surface in
    the all-transcript union -- this is exactly the 5227172 case."""
    vep_results = [make_vep_result(
        chrom="11", pos=5227172, ref="G", alt="A",
        transcript_consequences=[
            {"transcript_id": CANONICAL_TX, "gene_id": TARGET_GENE,
             "consequence_terms": ["upstream_gene_variant"]},
            {"transcript_id": "ENST00000000001", "gene_id": TARGET_GENE,
             "consequence_terms": ["5_prime_UTR_variant"]},
        ],
    )]

    df = parse_vep_results(vep_results, CANONICAL_TX, TARGET_GENE)
    row = df.iloc[0]

    assert row["vep_canonical_terms"] == ["upstream_gene_variant"], (
        "Canonical-only term should still be tracked separately for comparison"
    )
    assert "5_prime_UTR_variant" in row["vep_all_transcript_terms"], (
        "5_prime_UTR_variant from the non-canonical transcript was lost from the union -- "
        "this is the exact bug that mislabeled the real 5227172 variant"
    )


def test_canonical_and_union_disagreement_is_visible():
    """When canonical and union terms disagree, both must be queryable -- not silently
    resolved by picking one and discarding the other."""
    vep_results = [make_vep_result(
        chrom="11", pos=5227172, ref="G", alt="A",
        transcript_consequences=[
            {"transcript_id": CANONICAL_TX, "gene_id": TARGET_GENE,
             "consequence_terms": ["upstream_gene_variant"]},
            {"transcript_id": "ENST00000000001", "gene_id": TARGET_GENE,
             "consequence_terms": ["5_prime_UTR_variant"]},
        ],
    )]

    df = parse_vep_results(vep_results, CANONICAL_TX, TARGET_GENE)
    row = df.iloc[0]

    assert set(row["vep_canonical_terms"]) != set(row["vep_all_transcript_terms"]), (
        "This fixture is constructed to disagree -- if this assertion fails, "
        "the fixture itself is wrong, not the code"
    )
