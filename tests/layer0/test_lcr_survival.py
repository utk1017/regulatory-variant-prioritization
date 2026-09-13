"""Locks in Problem #1 from the review: a strict VEP-term-first filter would have dropped
deep-intronic/intergenic enhancers like the HBB LCR. Layer 0 must flag them via interval overlap
alone, independent of any VEP consequence term."""
from regvar.layer0.intervals import flag_interval_overlap
from regvar.layer0.pipeline import build_evidence


def test_lcr_like_variant_is_flagged_via_interval_overlap(synthetic_variants_df, synthetic_reg_intervals_df):
    result = flag_interval_overlap(synthetic_variants_df, synthetic_reg_intervals_df)

    lcr_row = result[result["pos"] == 5280000].iloc[0]
    assert "dELS" in lcr_row["interval_evidence"], (
        "Variant inside the LCR-like interval was not captured by flag_interval_overlap"
    )

    other_row = result[result["pos"] == 5225000].iloc[0]
    assert other_row["interval_evidence"] == [], (
        "Variant outside any interval should have empty interval_evidence, not a false positive"
    )


def test_lcr_like_variant_survives_to_is_regulatory(synthetic_variants_df, synthetic_reg_intervals_df):
    result = flag_interval_overlap(synthetic_variants_df, synthetic_reg_intervals_df)
    result["vep_all_transcript_terms"] = [[] for _ in range(len(result))]
    result["vep_regulatory_terms"] = [[] for _ in range(len(result))]

    evidence = result.apply(build_evidence, axis=1)
    is_regulatory = evidence.str.len() > 0

    lcr_idx = result.index[result["pos"] == 5280000][0]
    assert is_regulatory[lcr_idx], (
        "LCR-like variant must be regulatory via interval evidence alone, "
        "with zero VEP consequence terms -- this is the exact case v1's VEP-term-first "
        "filter would have dropped"
    )
