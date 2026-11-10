import json
import pandas as pd
import pytest
from regvar.layer1.variants import Variant, canon_key
from regvar.layer2.schema import validate_layer2_availability_output
from regvar.layer2.tiers import (
    compute_evidence_tiers,
    label_coverage_warnings,
    tier_balance_warnings,
)


def test_stratification_tier_boundaries():
    """Verify correct assignment of Tier 1 (>=2), Tier 2 (1), and Tier 3 (0)."""
    v_rich = Variant("11", 501, "A", "T", canon_key("11", 501, "A", "T"))
    v_mod = Variant("11", 502, "A", "T", canon_key("11", 502, "A", "T"))
    v_sparse = Variant("11", 503, "A", "T", canon_key("11", 503, "A", "T"))

    ev_by_key = {
        # Rich: eQTL + MPRA tested; GWAS is position-level and does not count toward allele-specific tier
        v_rich.key: [
            {
                "modality": "gwas", "missingness_flag": "observed_positive",
                "evidence_value": json.dumps({"matching_level": "position", "allele_specific": False}),
            },
            {"modality": "eqtl", "missingness_flag": "observed_null", "evidence_value": "{}"},
            {"modality": "mpra", "missingness_flag": "observed_positive", "evidence_value": "{}"},
        ],
        # Moderate: only eQTL tested
        v_mod.key: [
            {"modality": "gwas", "missingness_flag": "absent", "evidence_value": "{}"},
            {"modality": "eqtl", "missingness_flag": "observed_positive", "evidence_value": "{}"},
        ],
        # Sparse: all absent
        v_sparse.key: [
            {"modality": "gwas", "missingness_flag": "absent", "evidence_value": "{}"},
            {"modality": "eqtl", "missingness_flag": "absent", "evidence_value": "{}"},
        ],
    }
    df = compute_evidence_tiers([v_rich, v_mod, v_sparse], ev_by_key)
    tier_map = dict(zip(df["variant_key"], df["evidence_tier"]))

    key_rich = f"{v_rich.chrom}:{v_rich.pos}:{v_rich.ref}:{v_rich.alt}"
    key_mod = f"{v_mod.chrom}:{v_mod.pos}:{v_mod.ref}:{v_mod.alt}"
    key_sparse = f"{v_sparse.chrom}:{v_sparse.pos}:{v_sparse.ref}:{v_sparse.alt}"

    assert tier_map[key_rich] == "Tier 1 (Evidence-Rich)"
    assert tier_map[key_mod] == "Tier 2 (Moderate Evidence)"
    assert tier_map[key_sparse] == "Tier 3 (Evidence-Sparse)"


def test_tier_and_label_warnings():
    df = pd.DataFrame({
        "evidence_tier": ["Tier 3 (Evidence-Sparse)"] * 6,
        "label": ["pathogenic"] * 6,
    })
    tw = tier_balance_warnings(df)
    assert any(w.startswith("Tier 1") and "0 variants" in w for w in tw)
    assert any(w.startswith("Tier 2") and "0 variants" in w for w in tw)
    lw = label_coverage_warnings(df)
    assert any("no 'benign'" in w for w in lw)
    assert label_coverage_warnings(df.drop(columns=["label"]))[0].startswith("No ClinVar")


def test_labels_matching_no_variant_does_not_crash():
    v = Variant("11", 700, "A", "G", canon_key("11", 700, "A", "G"))
    df = compute_evidence_tiers(
        [v],
        {v.key: []},
        clinvar_labels={("11", 1, "A", "G"): {"CLNSIG": "x", "label": "benign", "review_stars": 1}},
    )
    assert "label" in df.columns and df["label"].isna().all()


def test_gnomad_af_status_is_validated():
    """validate_layer2_availability_output validates all status columns including gnomad_af_status."""
    v = Variant("11", 904, "C", "T", canon_key("11", 904, "C", "T"))
    ev_by_key = {
        v.key: [{
            "modality": "gnomad_af", "context": "gnomad_r4:genomes",
            "missingness_flag": "observed_positive", "evidence_value": json.dumps({"af": 0.2}),
        }]
    }
    df = compute_evidence_tiers([v], ev_by_key)
    validate_layer2_availability_output(df)
    df.loc[0, "gnomad_af_status"] = "not_a_real_status"
    with pytest.raises(ValueError, match="Invalid status value"):
        validate_layer2_availability_output(df)
