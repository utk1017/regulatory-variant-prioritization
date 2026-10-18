import json
import time
import pandas as pd
from regvar.layer1.gnomad import (
    build_gnomad_records,
    collect_gnomad_evidence,
    maf_from_gnomad,
)
from regvar.layer1.variants import to_variants


def _variants(rows):
    return to_variants(pd.DataFrame([{"chrom": "11", "pos": p, "ref": r, "alt": a} for p, r, a in rows]))


def test_gnomad_collection_and_af_calculation():
    variants = _variants([(100, "A", "G"), (200, "C", "T"), (900000, "G", "A")])

    def mock_fetch(method, url, params=None, json_body=None):
        variables = (json_body or {}).get("variables", {})
        start = variables.get("start", 0)
        if start > 800000:
            return None, "http_503", False
        return {
            "data": {
                "region": {
                    "variants": [
                        {
                            "variant_id": "11-100-A-G",
                            "pos": 100,
                            "ref": "A",
                            "alt": "G",
                            "genome": {"ac": 3, "an": 150000},
                        },
                        {
                            "variant_id": "11-150-A-G",
                            "pos": 150,
                            "ref": "A",
                            "alt": "G",
                            "genome": None,
                        },
                    ]
                }
            }
        }, None, False

    gnomad = collect_gnomad_evidence(variants, "gnomad_r4", window_bp=25000, delay_s=0.0, fetch=mock_fetch)
    assert len(gnomad["windows"]) == 2
    assert len(gnomad["by_key"]) == 1

    maf = maf_from_gnomad(gnomad)
    assert abs(maf[("11", 100, "A", "G")] - 3 / 150000) < 1e-12

    records = build_gnomad_records(gnomad, variants, "gnomad_r4")
    recs_by_pos = {r.pos: r for r in records}

    assert recs_by_pos[100].missingness_flag == "observed_positive"
    assert json.loads(recs_by_pos[100].evidence_value)["ac"] == 3

    assert recs_by_pos[200].missingness_flag == "absent"
    assert json.loads(recs_by_pos[200].evidence_value)["reason"] == "not_observed_in_genome_callset"

    assert recs_by_pos[900000].missingness_flag == "absent"
    assert json.loads(recs_by_pos[900000].evidence_value)["reason"] == "query_failed"


def test_gnomad_cache_hits_skip_sleep():
    """Verify that when from_cache is True, rate-limiting sleep is skipped."""
    variants = _variants([(100, "A", "G"), (200000, "C", "T")])

    def mock_cached_fetch(method, url, params=None, json_body=None):
        return {"data": {"region": {"variants": []}}}, None, True  # from_cache=True

    t0 = time.time()
    collect_gnomad_evidence(variants, "gnomad_r4", window_bp=25000, delay_s=5.0, fetch=mock_cached_fetch)
    elapsed = time.time() - t0
    assert elapsed < 0.5, f"Cache hits must skip delay sleep; took {elapsed:.2f}s"
