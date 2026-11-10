import json
import os
import tempfile
import pandas as pd
import pytest
import yaml

from regvar.layer2.pipeline import run_layer2
from regvar.layer2.schema import (
    LAYER2_AGREEMENT_COLUMNS,
    LAYER2_AVAILABILITY_COLUMNS,
    LAYER2_TRIGGER_COLUMNS,
)


def test_run_layer2_synthetic():
    with tempfile.TemporaryDirectory() as tmpdir:
        cfg_path = os.path.join(tmpdir, "test_config.yaml")
        config_data = {
            "gene_symbol": "HBB_TEST",
            "layer2": {
                "target_tissues": ["Whole_Blood"],
                "target_cell_types": ["K562"],
                "strict": False,
            },
        }
        with open(cfg_path, "w", encoding="utf-8") as fh:
            yaml.dump(config_data, fh)

        avail_df, trig_df, agree_df = run_layer2(cfg_path)

        assert not avail_df.empty
        assert not trig_df.empty
        assert not agree_df.empty

        for col in LAYER2_AVAILABILITY_COLUMNS:
            assert col in avail_df.columns
        for col in LAYER2_TRIGGER_COLUMNS:
            assert col in trig_df.columns
        for col in LAYER2_AGREEMENT_COLUMNS:
            assert col in agree_df.columns

        # Verify output manifest in synthetic_demo
        manifest_path = os.path.join("data", "outputs", "HBB_TEST", "synthetic_demo", "layer2_manifest.json")
        assert os.path.exists(manifest_path)
        with open(manifest_path, "r", encoding="utf-8") as fh:
            manifest = json.load(fh)
        assert manifest["is_synthetic"] is True
        assert manifest["layer"] == 2
        assert manifest["gene_symbol"] == "HBB_TEST"


def test_run_layer2_with_live_inputs():
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create synthetic Layer 1 CSV
        l1_csv = os.path.join(tmpdir, "mock_layer1_evidence.csv")
        records = [
            {
                "chrom": "11", "pos": 5227000, "ref": "C", "alt": "T",
                "modality": "eqtl", "source_type": "experimental",
                "context": "Whole_Blood", "missingness_flag": "observed_positive",
                "evidence_value": json.dumps({"beta": -0.4, "pvalue": 1e-5}),
            },
            {
                "chrom": "11", "pos": 5227000, "ref": "C", "alt": "T",
                "modality": "gwas", "source_type": "experimental",
                "context": None, "missingness_flag": "observed_null",
                "evidence_value": json.dumps({"reason": "subthreshold_association"}),
            },
            {
                "chrom": "11", "pos": 5227000, "ref": "C", "alt": "T",
                "modality": "mpra", "source_type": "experimental",
                "context": "K562", "missingness_flag": "observed_positive",
                "evidence_value": json.dumps({"log2_allelic_skew": 1.2, "pvalue": 1e-4}),
            },
            {
                "chrom": "11", "pos": 5227000, "ref": "C", "alt": "T",
                "modality": "crispr", "source_type": "experimental",
                "context": None, "missingness_flag": "absent",
                "evidence_value": json.dumps({"reason": "no_data_source_configured"}),
            },
            {
                "chrom": "11", "pos": 5227000, "ref": "C", "alt": "T",
                "modality": "gnomad_af", "source_type": "curated_clinical",
                "context": "genome", "missingness_flag": "observed_positive",
                "evidence_value": json.dumps({"af": 0.05}),
            },
        ]
        pd.DataFrame(records).to_csv(l1_csv, index=False)

        cfg_path = os.path.join(tmpdir, "test_config.yaml")
        config_data = {
            "gene_symbol": "HBB_MOCK",
            "layer1_csv_path": l1_csv,
            "layer2": {
                "target_tissues": ["Whole_Blood"],
                "target_cell_types": ["K562"],
                "strict": True,
            },
        }
        with open(cfg_path, "w", encoding="utf-8") as fh:
            yaml.dump(config_data, fh)

        avail_df, trig_df, agree_df = run_layer2(cfg_path)

        assert len(avail_df) == 1
        assert avail_df.iloc[0]["eqtl_status"] == "observed_positive"
        assert avail_df.iloc[0]["evidence_tier"] == "Tier 1 (Evidence-Rich)"

        # Triggers should cover crispr in K562
        assert len(trig_df) == 1
        assert set(trig_df["modality"]) == {"crispr"}
        assert set(trig_df["context"]) == {"K562"}

        # Agreement check should cover observed eQTL and observed MPRA
        assert len(agree_df) == 2
        assert set(agree_df["modality"]) == {"eqtl", "mpra"}

        # Verify live manifest
        manifest_path = os.path.join("data", "outputs", "HBB_MOCK", "layer2_manifest.json")
        assert os.path.exists(manifest_path)
        with open(manifest_path, "r", encoding="utf-8") as fh:
            manifest = json.load(fh)
        assert manifest["is_synthetic"] is False


def test_run_layer2_strict_missing_raises():
    with tempfile.TemporaryDirectory() as tmpdir:
        cfg_path = os.path.join(tmpdir, "test_config.yaml")
        config_data = {
            "gene_symbol": "NONEXISTENT_GENE_12345",
            "layer2": {
                "strict": True,
            },
        }
        with open(cfg_path, "w", encoding="utf-8") as fh:
            yaml.dump(config_data, fh)

        with pytest.raises(FileNotFoundError, match="Strict mode enabled"):
            run_layer2(cfg_path)
