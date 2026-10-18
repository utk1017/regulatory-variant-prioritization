import json
import os
import tempfile
import pandas as pd
import pytest

from regvar.layer1.local_tables import (
    CRISPR_COLUMNS,
    MPRA_COLUMNS,
    build_local_table_records,
    load_local_evidence_table,
    mpra_value,
    crispr_value,
)
from regvar.layer1.variants import to_variants


def _variants(rows):
    return to_variants(pd.DataFrame([{"chrom": "11", "pos": p, "ref": r, "alt": a} for p, r, a in rows]))


def test_local_table_classification():
    variants = _variants([(100, "A", "G"), (200, "C", "T"), (300, "G", "A"), (400, "T", "C")])

    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "mpra.csv")
        pd.DataFrame([
            {"chrom": "11", "pos": 100, "ref": "A", "alt": "G", "cell_type": "K562", "log2_allelic_skew": 1.8, "pvalue": 0.001, "dataset_name": "screen1"},
            {"chrom": "11", "pos": 200, "ref": "C", "alt": "T", "cell_type": "K562", "log2_allelic_skew": 0.1, "pvalue": 0.4, "dataset_name": "screen1"},
            {"chrom": "11", "pos": 300, "ref": "G", "alt": "A", "cell_type": "K562", "log2_allelic_skew": 0.1, "pvalue": None, "dataset_name": "screen1"},
        ]).to_csv(path, index=False)

        idx = load_local_evidence_table(path, MPRA_COLUMNS, "MPRA", ("log2_allelic_skew",))
        recs = build_local_table_records("mpra", variants, idx, p_threshold=0.05, value_fn=mpra_value)
        by_pos = {r.pos: r for r in recs}

        assert by_pos[100].missingness_flag == "observed_positive"
        assert by_pos[200].missingness_flag == "observed_null"
        assert by_pos[300].missingness_flag == "absent"
        assert json.loads(by_pos[300].evidence_value)["reason"] == "p_value_missing"
        assert by_pos[400].missingness_flag == "absent"
        assert json.loads(by_pos[400].evidence_value)["reason"] == "not_tested_in_provided_table"


def test_local_table_errors():
    with tempfile.TemporaryDirectory() as tmpdir:
        bad_path = os.path.join(tmpdir, "bad.csv")
        pd.DataFrame([{"chrom": "11", "pos": 100}]).to_csv(bad_path, index=False)

        with pytest.raises(ValueError, match="missing expected column"):
            load_local_evidence_table(bad_path, MPRA_COLUMNS, "MPRA")

        with pytest.raises(FileNotFoundError):
            load_local_evidence_table(os.path.join(tmpdir, "nonexistent.csv"), MPRA_COLUMNS, "MPRA")


def test_unconfigured_local_table_returns_absent():
    variants = _variants([(100, "A", "G")])
    recs = build_local_table_records("crispr", variants, index=None, p_threshold=0.05, value_fn=crispr_value)
    assert len(recs) == 1
    assert recs[0].missingness_flag == "absent"
    assert json.loads(recs[0].evidence_value)["reason"] == "no_data_source_configured"
