import contextlib
import json
import os
from pathlib import Path
import tempfile
from unittest import mock
import zipfile

import pandas as pd
import pytest
import yaml

from regvar.layer0.intervals import flag_interval_overlap
from regvar.layer0.vep import parse_vep_results
from regvar.layer1.cache import HTTPCache, make_fetch
from regvar.layer1.gwas import (
    build_gwas_records,
    collect_gwas_evidence,
    collect_gwas_evidence_from_download,
    collect_gwas_with_fallback,
    fetch_paginated_hal,
)
from regvar.layer1.intake import load_layer0_output
from regvar.layer1.schema import records_to_dataframe
from regvar.layer1.variants import Variant, canon_key
from regvar.layer2.pipeline import (
    audit_context_scope,
    enforce_data_quality_gate,
    scan_query_errors,
    tier1_feasibility,
)


@contextlib.contextmanager
def _in_tmpdir():
    old = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
        os.chdir(d)
        try:
            yield Path(d)
        finally:
            os.chdir(old)


class _FakeResp:
    def __init__(self, status, payload=None):
        self.status_code, self._p, self.headers = status, payload, {}

    def json(self):
        if self._p is None:
            raise ValueError("no json")
        return self._p


def ensure_project_config(path="configs/hbb.yaml"):
    p = Path(path)
    if p.exists():
        return False
    p.parent.mkdir(parents=True, exist_ok=True)
    default_yaml = """gene_symbol: HBB
flank_bp: 100000
layer0:
  clinvar_snapshot_date: "20261019"
layer1:
  gwas_source: auto
layer2:
  target_cell_types: [HEL92.1.7]
"""
    p.write_text(default_yaml, encoding="utf-8")
    return True


def test_default_config_is_created_once_and_parses():
    with _in_tmpdir():
        assert ensure_project_config("configs/t.yaml") is True
        cfg = yaml.safe_load(Path("configs/t.yaml").read_text())
        assert cfg["layer2"]["target_cell_types"] == ["HEL92.1.7"]
        assert cfg["layer0"]["clinvar_snapshot_date"] == "20261019"
        Path("configs/t.yaml").write_text("gene_symbol: EDITED\n")
        assert ensure_project_config("configs/t.yaml") is False
        assert "EDITED" in Path("configs/t.yaml").read_text()


def test_http_410_is_retried_then_reported():
    calls = {"n": 0}

    def always_410(*a, **k):
        calls["n"] += 1
        return _FakeResp(410)

    cache = HTTPCache(cache_dir=None, enabled=False)
    f = make_fetch(cache, max_retries=3, backoff=1.0)
    with mock.patch.object(f.session, "request", always_410), mock.patch("time.sleep"):
        data, err, _ = f("GET", "https://example.test/x")
    assert err == "http_410" and data is None and calls["n"] == 3

    seq = iter([_FakeResp(410), _FakeResp(200, {"ok": 1})])
    f2 = make_fetch(cache, max_retries=3, backoff=1.0)
    with mock.patch.object(f2.session, "request", lambda *a, **k: next(seq)), mock.patch("time.sleep"):
        data, err, _ = f2("GET", "https://example.test/y")
    assert err is None and data == {"ok": 1}


def test_gwas_response_shape_tolerance_and_guard():
    hal = {"_embedded": {"snps": [{"rs_id": "rs1"}]}, "_links": {"next": {"href": "https://x/next{?page}"}}}
    plain = {"content": [{"rs_id": "rs2"}], "page": {"number": 0}}
    empty_valid = {"_links": {}, "page": {"totalElements": 0}}
    weird = {"unexpected": "shape"}
    pages = {
        "https://x/start": hal,
        "https://x/next": {"_embedded": {"snps": [{"rs_id": "rs3"}]}, "_links": {}},
    }
    fetch = lambda m, url, params=None: (pages[url], None, False)
    items, err = fetch_paginated_hal("https://x/start", {}, "snps", fetch)
    assert err is None and [i["rs_id"] for i in items] == ["rs1", "rs3"]
    assert fetch_paginated_hal("https://x/p", {}, "snps", lambda m, u, params=None: (plain, None, False)) == ([{"rs_id": "rs2"}], None)
    assert fetch_paginated_hal("https://x/e", {}, "snps", lambda m, u, params=None: (empty_valid, None, False)) == ([], None)
    items, err = fetch_paginated_hal("https://x/w", {}, "snps", lambda m, u, params=None: (weird, None, False))
    assert err == "unrecognized_response_shape" and items == []


def test_gwas_failure_propagates_to_layer2_and_gate():
    v = Variant("11", 5225600, "A", "G", canon_key("11", 5225600, "A", "G"))
    failing = lambda m, url, params=None, json_body=None: (None, "http_410", False)
    gwas = collect_gwas_evidence("11", [v.pos], fetch=failing, window_bp=25000)
    recs = build_gwas_records(gwas, [v], 5e-8)
    assert recs[0].missingness_flag == "absent" and "incomplete_region_query" in recs[0].evidence_value
    df = records_to_dataframe(recs)
    warns = scan_query_errors(df)
    assert len(warns) == 1 and "gwas" in warns[0] and "http_410" in warns[0]
    with pytest.raises(RuntimeError, match="Layer 2 halted"):
        enforce_data_quality_gate(warns, halt=True)
    assert enforce_data_quality_gate(warns, halt=False) == "DEGRADED"
    assert enforce_data_quality_gate([], halt=True) == "OK"


def test_layer0_manifest_is_discovered_and_gene_meta_pinned():
    with _in_tmpdir():
        Path("layer0_outputs").mkdir()
        pd.DataFrame({"chrom": ["11"], "pos": [5227002], "ref": ["T"], "alt": ["C"], "label": ["benign"]}).to_csv(
            "layer0_outputs/HBB_layer0_regulatory_truthset.csv", index=False
        )
        gm = {
            "gene_symbol": "HBB",
            "gene_id": "ENSG00000244734",
            "chrom": "11",
            "start": 5225464,
            "end": 5229395,
            "strand": -1,
        }
        Path("layer0_outputs/HBB_layer0_manifest.json").write_text(
            json.dumps({"gene_metadata": gm, "clinvar_snapshot_date": "20261019"})
        )
        df, meta, info = load_layer0_output(
            "HBB", None, fetch=lambda *a, **k: (_ for _ in ()).throw(AssertionError("no network allowed"))
        )
        assert meta == gm and info["manifest_path"] and info["clinvar_snapshot_date"] == "20261019"


def test_scope_audit_catches_wrong_cell_type():
    ev_df = pd.DataFrame([
        {
            "chrom": "11",
            "pos": 1,
            "ref": "A",
            "alt": "G",
            "modality": "mpra",
            "source_type": "experimental",
            "evidence_value": "{}",
            "context": "K562",
            "missingness_flag": "observed_positive",
        },
        {
            "chrom": "11",
            "pos": 2,
            "ref": "A",
            "alt": "G",
            "modality": "mpra",
            "source_type": "experimental",
            "evidence_value": "{}",
            "context": None,
            "missingness_flag": "absent",
        },
    ])
    warns, tbl = audit_context_scope(ev_df, ["Whole_Blood"], ["HEL92.1.7"])
    assert any("NONE of the mpra" in w for w in warns) and int(tbl.loc[tbl.modality == "mpra", "records_excluded"].iloc[0]) == 1
    warns_ok, _ = audit_context_scope(ev_df, ["Whole_Blood"], ["K562"])
    assert not warns_ok


def test_tier1_feasibility_lists_promotable_variants():
    av = pd.DataFrame({
        "variant_key": ["a", "b", "c", "d"],
        "eqtl_status": ["absent", "observed_positive", "absent", "absent"],
        "mpra_status": ["observed_positive", "observed_positive", "absent", "observed_null"],
        "crispr_status": ["absent", "absent", "absent", "absent"],
        "evidence_tier": [
            "Tier 2 (Moderate Evidence)",
            "Tier 1 (Evidence-Rich)",
            "Tier 3 (Evidence-Sparse)",
            "Tier 2 (Moderate Evidence)",
        ],
    })
    c = tier1_feasibility(av, 2)
    assert list(c["variant_key"]) == ["a", "d"] and c.iloc[0]["has_modalities"] == "mpra"


def test_layer0_logic_offline():
    vep = [
        {
            "start": 5226000,
            "seq_region_name": "11",
            "allele_string": "T/C",
            "transcript_consequences": [
                {
                    "gene_id": "ENSG00000244734",
                    "transcript_id": "ENST00000335295",
                    "consequence_terms": ["5_prime_UTR_variant"],
                },
                {"gene_id": "OTHER", "consequence_terms": ["missense_variant"]},
            ],
        },
        {
            "start": 5226100,
            "seq_region_name": "11",
            "allele_string": "A/G",
            "transcript_consequences": [],
            "regulatory_feature_consequences": [],
        },
    ]
    vdf = parse_vep_results(vep, "ENST00000335295", "ENSG00000244734")
    assert list(vdf["vep_all_transcript_terms"][0]) == ["5_prime_UTR_variant"]

    base = pd.DataFrame({
        "chrom": ["11", "11"],
        "pos": [5226000, 5226100],
        "ref": ["T", "A"],
        "alt": ["C", "G"],
        "CLNSIG": ["Pathogenic", "Benign"],
        "CLNREVSTAT": ["x", "y"],
        "review_stars": [2, 1],
        "label": ["pathogenic", "benign"],
        "interval_evidence": [[], []],
    })
    flagged = flag_interval_overlap(
        base[["chrom", "pos", "ref", "alt"]].assign(label=["pathogenic", "benign"]),
        pd.DataFrame({"Chromosome": ["11"], "Start": [5225990], "End": [5226010], "ccre_label": ["promoter-like"]}),
    )
    assert flagged["interval_evidence"].tolist() == [["promoter-like"], []]


def _make_bulk_zip(path):
    header = [
        "DATE ADDED TO CATALOG", "PUBMEDID", "DISEASE/TRAIT", "CHR_ID", "CHR_POS",
        "MAPPED_GENE", "SNPS", "P-VALUE", "PVALUE_MLOG", "MAPPED_TRAIT", "STUDY ACCESSION"
    ]
    rows = [
        ["2020-01-01", "1", "Fetal hemoglobin", "11", "5227002", "HBB", "rs1", "5E-8", "7.3", "fetal hemoglobin measurement", "GCST1"],
        ["2020-01-01", "2", "Hematocrit", "11; 11", "5227002; 5227100", "HBB", "rs1; rs2", "", "9", "hematocrit", "GCST2"],
        ["2020-01-01", "3", "Interaction", "11 x 12", "5227500 x 100", "HBB", "rs9 x rs8", "1E-30", "29", "x", "GCST3"],
        ["2020-01-01", "4", "Other chromosome", "12", "5227002", "X", "rs7", "1E-30", "29", "y", "GCST4"],
        ["2020-01-01", "5", "Weak", "11", "5300000", "HBB", "rs5", "3E-3", "2.5", "weak trait", "GCST5"],
    ]
    tsv = "\t".join(header) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("gwas-catalog-associations_ontology-annotated.tsv", tsv)


def test_gwas_bulk_file_parsing_and_records():
    with _in_tmpdir() as d:
        zp = str(d / "bulk.zip")
        _make_bulk_zip(zp)
        positions = [5227002, 5227100, 5228000, 5300000]
        res = collect_gwas_evidence_from_download("11", positions, zp, window_bp=25000)
        assert all(w["error"] is None for w in res["windows"]) and res["n_snps_in_region"] == 4
        assert set(res["by_pos"]) == {5227002, 5227100, 5300000}
        e = res["by_pos"][5227002]
        assert e["rs_ids"] == ["rs1"] and len(e["assocs"]) == 2
        pv = sorted(a["pvalue"] for a in e["assocs"])
        assert abs(pv[0] - 1e-9) < 1e-12 and abs(pv[1] - 5e-8) < 1e-12

        vs = [Variant("11", p, "A", "G", canon_key("11", p, "A", "G")) for p in positions]
        recs = {r.pos: r for r in build_gwas_records(res, vs, 5e-8)}
        assert recs[5227002].missingness_flag == "observed_positive" and recs[5227100].missingness_flag == "observed_positive"
        assert recs[5300000].missingness_flag == "observed_null"
        assert json.loads(recs[5228000].evidence_value)["reason"] == "not_in_catalog"
        assert "bulk download" in json.loads(recs[5227002].evidence_value)["source"]

        bad = str(d / "bad.zip")
        Path(bad).write_text("not a zip")
        broken = collect_gwas_evidence_from_download("11", positions, bad)
        assert broken["windows"][0]["error"].startswith("bulk_file_unreadable")


def test_gwas_auto_fallback_and_modes():
    pos = [5227002]
    failing = lambda m, u, params=None, json_body=None: (None, "http_410", False)
    good_bulk = {
        "windows": [{"chrom": "11", "start": 5227002, "stop": 5227002, "error": None, "n_snps": 1}],
        "by_pos": {5227002: {"rs_ids": ["rs1"], "assocs": [{"traits": ["t"], "pvalue": 1e-9, "rs_id": "rs1"}], "error": None}},
        "n_snps_in_region": 1,
        "source_label": "bulk",
    }
    probes = []
    r = collect_gwas_with_fallback(
        "11",
        pos,
        failing,
        {"gwas_source": "auto"},
        bulk_fn=lambda c, p, L: (good_bulk, {"sha256": "x"}, None),
        probe_fn=lambda *a: probes.append(a),
    )
    assert r["used"] == "download" and not r["warnings"] and r["notes"] and probes and r["file_info"] == {"sha256": "x"}

    r = collect_gwas_with_fallback(
        "11",
        pos,
        failing,
        {"gwas_source": "auto"},
        bulk_fn=lambda c, p, L: (None, None, "bulk_download_failed:X"),
        probe_fn=lambda *a: None,
    )
    assert r["used"] == "api" and any("also unavailable" in w for w in r["warnings"]) and any("REST API" in w for w in r["warnings"])

    called = []
    r = collect_gwas_with_fallback(
        "11", pos, failing, {"gwas_source": "api"}, bulk_fn=lambda *a: called.append(1), probe_fn=lambda *a: None
    )
    assert not called and r["warnings"] and r["used"] == "api"

    ok_fetch = lambda m, u, params=None, json_body=None: ({"_embedded": {"snps": []}}, None, False)
    r = collect_gwas_with_fallback(
        "11", pos, ok_fetch, {"gwas_source": "auto"}, bulk_fn=lambda *a: called.append(2), probe_fn=lambda *a: None
    )
    assert r["used"] == "api" and not r["warnings"] and called == []

    r = collect_gwas_with_fallback(
        "11", pos, failing, {"gwas_source": "download"}, bulk_fn=lambda c, p, L: (None, None, "bulk_download_failed:X")
    )
    assert r["used"] == "download" and r["warnings"]

    with pytest.raises(ValueError, match="layer1.gwas_source"):
        collect_gwas_with_fallback("11", pos, failing, {"gwas_source": "nope"})
