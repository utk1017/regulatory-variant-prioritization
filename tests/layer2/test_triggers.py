import json
import pandas as pd
import pytest
from regvar.layer1.variants import Variant, canon_key
from regvar.layer2.schema import (
    AG_OUTPUT_TYPES,
    MODALITY_META,
    validate_layer2_triggers_output,
)
from regvar.layer2.triggers import (
    compile_agreement_check_queue,
    compile_alphagenome_triggers,
)


def test_transient_error_handling():
    """Query failures must NOT trigger AlphaGenome and must be flagged as warnings."""
    v = Variant("11", 100, "A", "G", canon_key("11", 100, "A", "G"))
    ev_by_key = {
        v.key: [{
            "modality": "eqtl",
            "context": "Whole_Blood",
            "missingness_flag": "absent",
            "evidence_value": json.dumps({"reason": "query_failed", "error": "HTTP 500 Server Error"}),
        }]
    }
    trigs, warns = compile_alphagenome_triggers(
        variants=[v],
        evidence_by_key=ev_by_key,
        completable_modalities=["eqtl"],
        target_tissues=["Whole_Blood"],
        head_mapping=MODALITY_META,
        gene_meta={"gene_symbol": "TEST"},
    )
    assert len(trigs) == 0, "Transient query failure must not generate an AlphaGenome trigger!"
    assert len(warns) == 1, "Transient query failure must be recorded in warnings!"


def test_tested_nulls_never_trigger_alphagenome():
    """Observed negative results (null) must never be filled by AlphaGenome."""
    v = Variant("11", 200, "C", "T", canon_key("11", 200, "C", "T"))
    ev_by_key = {
        v.key: [{
            "modality": "eqtl",
            "context": "Whole_Blood",
            "missingness_flag": "observed_null",
            "evidence_value": json.dumps({"reason": "common_variant_not_significant", "maf": 0.25}),
        }]
    }
    trigs, _ = compile_alphagenome_triggers(
        variants=[v],
        evidence_by_key=ev_by_key,
        completable_modalities=["eqtl"],
        target_tissues=["Whole_Blood"],
        head_mapping=MODALITY_META,
        gene_meta={"gene_symbol": "TEST"},
    )
    assert len(trigs) == 0, "Tested null must NEVER trigger AlphaGenome!"


def test_non_completable_modalities():
    """GWAS and gnomAD AF must NEVER be trigger targets for AlphaGenome."""
    v = Variant("11", 300, "G", "C", canon_key("11", 300, "G", "C"))
    ev_by_key = {
        v.key: [
            {"modality": "gwas", "missingness_flag": "absent", "evidence_value": json.dumps({"reason": "not_in_catalog"})},
            {"modality": "gnomad_af", "missingness_flag": "absent", "evidence_value": json.dumps({"reason": "not_found"})},
        ]
    }
    trigs, _ = compile_alphagenome_triggers(
        variants=[v],
        evidence_by_key=ev_by_key,
        completable_modalities=["eqtl", "mpra", "crispr"],
        target_tissues=["Whole_Blood"],
        head_mapping=MODALITY_META,
        gene_meta={"gene_symbol": "TEST"},
    )
    modalities_in_trigs = {t["modality"] for t in trigs}
    assert "gwas" not in modalities_in_trigs, "GWAS cannot be completed by AlphaGenome!"
    assert "gnomad_af" not in modalities_in_trigs, "gnomAD AF cannot be completed by AlphaGenome!"


def test_agreement_queue_strict_isolation():
    """Agreement checks must only contain observed measurements and never overlap with triggers."""
    v = Variant("11", 400, "T", "G", canon_key("11", 400, "T", "G"))
    ev_by_key = {
        v.key: [
            {"modality": "mpra", "context": "K562", "missingness_flag": "observed_positive", "evidence_value": json.dumps({"skew": 1.2})},
            {"modality": "crispr", "context": "K562", "missingness_flag": "absent", "evidence_value": json.dumps({"reason": "no_data_source_configured"})},
        ]
    }
    trigs, _ = compile_alphagenome_triggers(
        [v], ev_by_key, ["mpra", "crispr"], [], MODALITY_META, {"gene_symbol": "TEST"}, target_cell_types=["K562"]
    )
    agrees = compile_agreement_check_queue(
        [v], ev_by_key, ["mpra", "crispr"], target_cell_types=["K562"]
    )

    trig_mods = {t["modality"] for t in trigs}
    agree_mods = {a["modality"] for a in agrees}

    assert "mpra" in agree_mods and "mpra" not in trig_mods
    assert "crispr" in trig_mods and "crispr" not in agree_mods


def test_target_tissues_respect_layer1_scope():
    """Un-queried target tissue raises exactly one CONFIG warning and queues nothing."""
    v = Variant("11", 600, "A", "G", canon_key("11", 600, "A", "G"))
    ev_by_key = {v.key: [{
        "modality": "eqtl", "context": "Whole_Blood", "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "not_testable_rare_or_unknown_frequency"}),
    }]}
    trigs, warns = compile_alphagenome_triggers(
        [v], ev_by_key, ["eqtl"], ["Whole_Blood", "Liver"], MODALITY_META, {"gene_symbol": "TEST"}
    )
    assert {t["context"] for t in trigs} == {"Whole_Blood"}
    liver_warns = [w for w in warns if w.startswith("CONFIG:") and "Liver" in w]
    assert len(liver_warns) == 1


def test_mpra_crispr_triggers_use_target_cell_types():
    v = Variant("11", 602, "C", "T", canon_key("11", 602, "C", "T"))
    ev_by_key = {v.key: [{
        "modality": "mpra", "context": None, "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "no_data_source_configured"}),
    }]}
    trigs, _ = compile_alphagenome_triggers(
        [v], ev_by_key, ["mpra", "crispr"], [], MODALITY_META, {"gene_symbol": "TEST"}, target_cell_types=["K562"]
    )
    assert {t["context"] for t in trigs} == {"K562"}
    assert {t["modality"] for t in trigs} == {"mpra", "crispr"}


def test_trigger_output_types_and_serialization():
    v = Variant("11", 603, "G", "A", canon_key("11", 603, "G", "A"))
    ev_by_key = {v.key: [
        {"modality": "eqtl", "context": "Whole_Blood", "missingness_flag": "absent",
         "evidence_value": json.dumps({"reason": "not_testable_rare_or_unknown_frequency", "maf": None})},
        {"modality": "mpra", "context": None, "missingness_flag": "absent",
         "evidence_value": json.dumps({"reason": "no_data_source_configured"})},
    ]}
    trigs, _ = compile_alphagenome_triggers(
        [v], ev_by_key, ["eqtl", "mpra", "crispr"], ["Whole_Blood"], MODALITY_META, {"gene_symbol": "TEST"}
    )
    df = pd.DataFrame(trigs)
    validate_layer2_triggers_output(df)
    assert "target_head" not in df.columns
    assert set(df["alphagenome_binding_status"]) == {"unverified"}
    for s in df["alphagenome_output_types"]:
        assert set(s.split(";")) <= set(AG_OUTPUT_TYPES)
    for s in df["source_evidence_value"]:
        assert isinstance(json.loads(s), dict), "source_evidence_value must not be double JSON-encoded"
    assert set(df.loc[df.modality == "mpra", "alphagenome_proxy_quality"]) == {"indirect_proxy"}


def test_mpra_crispr_per_cell_type_granularity():
    """A hit in ONE target cell type must not suppress completion for a second never-tested cell type."""
    v = Variant("11", 900, "A", "G", canon_key("11", 900, "A", "G"))
    ev_by_key = {v.key: [
        {"modality": "mpra", "context": "K562", "missingness_flag": "observed_positive",
         "evidence_value": json.dumps({"log2_allelic_skew": 2.1, "pvalue": 1e-7})},
    ]}
    trigs, _ = compile_alphagenome_triggers(
        [v], ev_by_key, ["mpra"], [], MODALITY_META, {"gene_symbol": "TEST"}, target_cell_types=["K562", "HUDEP2"]
    )
    contexts = {t["context"] for t in trigs}
    assert "K562" not in contexts, "K562 has a real hit -- must not be re-triggered"
    assert "HUDEP2" in contexts, "HUDEP2 was never tested at all -- MUST be triggered"


def test_mpra_crispr_default_context_reason_preserved():
    """A context-less ('default') absent record must carry its real reason into triggers."""
    v = Variant("11", 901, "C", "T", canon_key("11", 901, "C", "T"))
    ev_by_key = {v.key: [
        {"modality": "mpra", "context": None, "missingness_flag": "absent",
         "evidence_value": json.dumps({"reason": "no_data_source_configured"})},
    ]}
    trigs, _ = compile_alphagenome_triggers(
        [v], ev_by_key, ["mpra"], [], MODALITY_META, {"gene_symbol": "TEST"}, target_cell_types=["K562", "HUDEP2"]
    )
    assert {t["context"] for t in trigs} == {"K562", "HUDEP2"}
    assert all(t["trigger_reason"] == "no_data_source_configured" for t in trigs)


def test_mpra_crispr_never_tested_cell_type_is_not_a_config_mismatch():
    """A target cell type absent entirely from the evidence is a gap, not a config mismatch."""
    v = Variant("11", 902, "G", "A", canon_key("11", 902, "G", "A"))
    ev_by_key = {v.key: [
        {"modality": "mpra", "context": "K562", "missingness_flag": "observed_null",
         "evidence_value": json.dumps({"reason": "tested_not_significant", "pvalue": 0.8})},
    ]}
    trigs, warns = compile_alphagenome_triggers(
        [v], ev_by_key, ["mpra"], [], MODALITY_META, {"gene_symbol": "TEST"}, target_cell_types=["K562", "HUDEP2"]
    )
    assert any(t["context"] == "HUDEP2" for t in trigs)
    assert not any("CONFIG" in w and "HUDEP2" in w for w in warns)
