import json
import os
import tempfile
import pandas as pd
import pytest
import yaml

from regvar.layer1.pipeline import run_layer1
from regvar.layer1.schema import (
    LAYER1_COLUMNS,
    make_record,
    records_to_dataframe,
    validate_layer1_output,
    validate_matrix_complete,
)
from regvar.layer1.variants import to_variants


def test_schema_validators():
    variants = to_variants(pd.DataFrame([{"chrom": "11", "pos": 100, "ref": "A", "alt": "G"}]))
    records = [
        make_record(variants[0], mod, "experimental", "absent", reason="test")
        for mod in ["gwas", "eqtl", "gnomad_af", "mpra", "crispr"]
    ]
    df = records_to_dataframe(records)
    validate_layer1_output(df)
    validate_matrix_complete(df, variants)

    # Incomplete matrix raises ValueError
    with pytest.raises(ValueError, match="Evidence matrix incomplete"):
        validate_matrix_complete(df.iloc[:2], variants)

    # Invalid column values raise ValueError
    bad_df = df.copy()
    bad_df.loc[0, "missingness_flag"] = "invalid_flag"
    with pytest.raises(ValueError, match="invalid missingness_flag"):
        validate_layer1_output(bad_df)


def test_run_layer1_end_to_end():
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create synthetic Layer 0 output
        gene_dir = os.path.join(tmpdir, "outputs", "HBB")
        os.makedirs(gene_dir, exist_ok=True)
        l0_csv = os.path.join(gene_dir, "layer0_regulatory_truthset.csv")
        pd.DataFrame([
            {"chrom": "11", "pos": 5227000, "ref": "C", "alt": "T"},
            {"chrom": "11", "pos": 5227100, "ref": "G", "alt": "A"},
        ]).to_csv(l0_csv, index=False)

        manifest_path = os.path.join(gene_dir, "manifest.json")
        with open(manifest_path, "w", encoding="utf-8") as fh:
            json.dump({
                "gene_metadata": {
                    "gene_symbol": "HBB",
                    "gene_id": "ENSG00000244734",
                    "chrom": "11",
                    "start": 5225464,
                    "end": 5229395,
                    "strand": -1,
                },
                "clinvar_snapshot_date": "2026-01-01",
            }, fh)

        # Create synthetic config
        cfg_path = os.path.join(tmpdir, "test_config.yaml")
        config_data = {
            "gene_symbol": "HBB",
            "layer1": {
                "layer0_csv_path": l0_csv,
                "http_cache_dir": os.path.join(tmpdir, "cache"),
                "use_http_cache": True,
                "gnomad_request_delay_s": 0.0,
            },
        }
        with open(cfg_path, "w", encoding="utf-8") as fh:
            yaml.dump(config_data, fh)

        # Mock fetch returning clean responses
        def mock_fetch(method, url, params=None, json_body=None):
            if "gnomad" in url:
                return {
                    "data": {
                        "region": {
                            "variants": [
                                {
                                    "variant_id": "11-5227000-C-T",
                                    "pos": 5227000,
                                    "ref": "C",
                                    "alt": "T",
                                    "genome": {"ac": 10, "an": 100000},
                                }
                            ]
                        }
                    }
                }, None, False
            if "singleNucleotidePolymorphisms" in url:
                return {"_embedded": {"singleNucleotidePolymorphisms": []}}, None, False
            if "reference/gene" in url:
                return {"data": [{"geneSymbol": "HBB", "gencodeId": "ENSG00000244734.11"}]}, None, False
            if "singleTissueEqtl" in url:
                return {"data": []}, None, False
            return None, "http_404", False

        # Run pipeline
        df = run_layer1(cfg_path, fetch=mock_fetch)
        assert df is not None
        # 2 variants x 5 modalities = 10 records
        assert len(df) == 10

        # Check output files were generated
        out_csv = os.path.join("data", "outputs", "HBB", "layer1_evidence.csv")
        out_manifest = os.path.join("data", "outputs", "HBB", "layer1_manifest.json")
        assert os.path.exists(out_csv)
        assert os.path.exists(out_manifest)

        with open(out_manifest, "r", encoding="utf-8") as fh:
            m = json.load(fh)
        assert m["layer"] == 1
        assert m["n_candidate_variants"] == 2
        assert m["n_evidence_records"] == 10
        assert "cache" in m
