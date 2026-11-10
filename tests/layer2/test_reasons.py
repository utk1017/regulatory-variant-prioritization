import json
import pytest
from regvar.layer1.variants import Variant, canon_key
from regvar.layer2.reasons import (
    _norm_ctx,
    aggregate_modality_status,
    classify_evidence_reason,
)
from regvar.layer2.tiers import compute_evidence_tiers
from regvar.layer2.triggers import compile_agreement_check_queue, compile_alphagenome_triggers
from regvar.layer2.schema import MODALITY_META


def test_multi_context_aggregation_order():
    """Observed positive in ANY context must override null or absent."""
    recs = [
        {"missingness_flag": "observed_null", "context": "Whole_Blood", "evidence_value": "{}"},
        {"missingness_flag": "observed_positive", "context": "Liver", "evidence_value": "{}"},
        {"missingness_flag": "absent", "context": "Brain", "evidence_value": "{}"},
    ]
    status, cat, _, _ = aggregate_modality_status(recs)
    assert status == "observed_positive", f"Expected positive, got {status}"
    assert cat == "POSITIVE"

    # Null must override absent
    recs_null = [
        {"missingness_flag": "observed_null", "context": "Whole_Blood", "evidence_value": "{}"},
        {"missingness_flag": "absent", "context": "Liver", "evidence_value": "{}"},
    ]
    status_null, cat_null, _, _ = aggregate_modality_status(recs_null)
    assert status_null == "observed_null", f"Expected null, got {status_null}"
    assert cat_null == "NULL_EXPERIMENTAL"


def test_nan_context_is_normalized():
    recs = [{
        "missingness_flag": "absent", "context": float("nan"),
        "evidence_value": json.dumps({"reason": "no_data_source_configured"}),
    }]
    _, _, _, by_ctx = aggregate_modality_status(recs)
    assert list(by_ctx) == ["default"], f"NaN context must normalize to 'default', got {list(by_ctx)}"
    assert _norm_ctx(None) == "default" and _norm_ctx("") == "default" and _norm_ctx("K562") == "K562"


def test_non_target_context_does_not_lift_tier():
    """A Liver eQTL hit must not make an erythroid gene's variant look evidence-rich."""
    v = Variant("11", 601, "A", "G", canon_key("11", 601, "A", "G"))
    ev_by_key = {v.key: [
        {"modality": "eqtl", "context": "Liver", "missingness_flag": "observed_positive", "evidence_value": "{}"},
        {"modality": "mpra", "context": "HepG2", "missingness_flag": "observed_positive", "evidence_value": "{}"},
        {"modality": "gwas", "context": None, "missingness_flag": "absent", "evidence_value": "{}"},
    ]}
    row = compute_evidence_tiers(
        [v],
        ev_by_key,
        regulatory_modalities=("gwas", "eqtl", "mpra", "crispr"),
        target_tissues=["Whole_Blood"],
        target_cell_types=["K562"],
    ).iloc[0]
    assert row["eqtl_status"] == "absent" and row["mpra_status"] == "absent"
    assert row["evidence_tier"] == "Tier 3 (Evidence-Sparse)"
    agrees = compile_agreement_check_queue(
        [v], ev_by_key, ["eqtl", "mpra"], target_tissues=["Whole_Blood"], target_cell_types=["K562"]
    )
    assert agrees == [], "Out-of-scope contexts must not enter the agreement queue"


def test_unknown_reason_is_never_completed_and_warns():
    """A reason string Layer 1 emits that Layer 2 doesn't recognize must not trigger AlphaGenome."""
    assert classify_evidence_reason("absent", reason="some_future_reason_layer1_added") == "UNKNOWN_REASON"
    assert classify_evidence_reason("absent", reason=None) == "GAP_BIOLOGICAL"

    v = Variant("11", 903, "A", "G", canon_key("11", 903, "A", "G"))
    ev_by_key = {v.key: [{
        "modality": "eqtl", "context": "Whole_Blood", "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "totally_new_reason_string"}),
    }]}
    trigs, warns = compile_alphagenome_triggers(
        [v], ev_by_key, ["eqtl"], ["Whole_Blood"], MODALITY_META, {"gene_symbol": "TEST"}
    )
    assert len(trigs) == 0, "an unrecognized reason must NOT silently generate an AlphaGenome trigger"
    assert any("unrecognized reason" in w for w in warns), "vocabulary drift must be reported in warnings"
