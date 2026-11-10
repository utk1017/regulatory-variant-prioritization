import json
import os
import tempfile
import pandas as pd
import pytest
from regvar.layer1.variants import canon_key, to_variants
from regvar.layer2.intake import (
    derive_fallback_gene_meta,
    generate_synthetic_hbb_evidence,
)
from regvar.layer2.pipeline import resolve_out_dir, write_layer2_manifest
from regvar.layer2.schema import validate_layer1_input


def test_empty_layer1_evidence_raises_clear_error():
    empty_df = pd.DataFrame(columns=[
        "chrom", "pos", "ref", "alt", "modality", "source_type",
        "evidence_value", "context", "missingness_flag",
    ])
    empty_variants = to_variants(empty_df)
    with pytest.raises(ValueError, match="zero candidate variants"):
        derive_fallback_gene_meta(empty_variants, "HBB")


def test_canon_key_matches_layer1_suffix_trimming():
    assert canon_key("11", 1000, "CAT", "GAT") == canon_key("11", 1000, "C", "G")
    assert canon_key("11", 5227100, "G", "GA") == canon_key("11", 5227099, "AG", "AGA")


def test_layer1_input_validation_rejects_bad_flags():
    bad_df = pd.DataFrame([{
        "chrom": "11", "pos": 1, "ref": "A", "alt": "G", "modality": "eqtl",
        "source_type": "experimental", "evidence_value": "{}", "context": None,
        "missingness_flag": "observed_postive",  # typo
    }])
    with pytest.raises(ValueError, match="invalid missingness_flag"):
        validate_layer1_input(bad_df)

    good_df = bad_df.copy()
    good_df.loc[0, "missingness_flag"] = "observed_positive"
    validate_layer1_input(good_df)


def test_synthetic_outputs_are_segregated():
    assert resolve_out_dir("HBB", False) == os.path.join("data", "outputs", "HBB")
    syn = resolve_out_dir("HBB", True)
    assert syn == os.path.join("data", "outputs", "HBB", "synthetic_demo")


def test_manifest_provenance_recording():
    with tempfile.TemporaryDirectory() as tmpdir:
        payload = {
            "gene_symbol": "HBB",
            "tier_distribution": {"Tier 1": 10, "Tier 3": 50},
            "alphagenome_triggers": {"total_jobs": 120},
        }
        path = write_layer2_manifest(tmpdir, payload)
        assert os.path.exists(path)
        with open(path, "r", encoding="utf-8") as f:
            loaded = json.load(f)
        assert loaded["gene_symbol"] == "HBB"
        assert "generated_at_utc" in loaded
        assert loaded["tier_distribution"]["Tier 3"] == 50
