import json
import pandas as pd
from regvar.layer1.gwas import build_gwas_records, collect_gwas_evidence
from regvar.layer1.variants import to_variants


def _variants(rows):
    return to_variants(pd.DataFrame([{"chrom": "11", "pos": p, "ref": r, "alt": a} for p, r, a in rows]))


def test_gwas_multi_window_failure_isolation():
    """Lock in Problem #2 from the review: a failure in one window must not relabel healthy windows."""
    variants = _variants([(100, "A", "G"), (50000, "C", "T"), (100000, "G", "A")])

    snp = lambda rs, pos: {
        "rsId": rs,
        "locations": [{"chromosomeName": "11", "chromosomePosition": pos}],
        "_links": {"snpAssociations": {"href": f"https://x/assoc/{rs}"}},
    }

    def mock_fetch(method, url, params=None, json_body=None):
        if "findByChromBpLocationRange" in url:
            start = params["bpStart"]
            if start < 1000:       # Window 1: fails
                return None, "http_500", False
            elif start < 60000:    # Window 2: healthy, empty
                return {"_embedded": {"singleNucleotidePolymorphisms": []}}, None, False
            else:                  # Window 3: healthy, has hit
                return {"_embedded": {"singleNucleotidePolymorphisms": [snp("rs3", 100000)]}}, None, False
        if "/assoc/" in url:
            return {"_embedded": {"associations": [{"pvalue": 1e-9, "efoTraits": [{"trait": "HbF"}]}]}}, None, False
        return None, "http_404", False

    gwas = collect_gwas_evidence("11", [100, 50000, 100000], fetch=mock_fetch, window_bp=25000)
    recs = build_gwas_records(gwas, variants, 5e-8)

    r100 = [r for r in recs if r.pos == 100][0]
    r50k = [r for r in recs if r.pos == 50000][0]
    r100k = [r for r in recs if r.pos == 100000][0]

    # Window 1 failed -> incomplete_region_query
    assert r100.missingness_flag == "absent"
    assert json.loads(r100.evidence_value)["reason"] == "incomplete_region_query"

    # Window 2 succeeded with no hits -> not_in_catalog (NOT incomplete_region_query!)
    assert r50k.missingness_flag == "absent"
    assert json.loads(r50k.evidence_value)["reason"] == "not_in_catalog"

    # Window 3 succeeded with significant hit -> observed_positive
    assert r100k.missingness_flag == "observed_positive"


def test_gwas_pagination_and_classification():
    variants = _variants([(100, "A", "G"), (300, "C", "T"), (500, "G", "A")])

    snp = lambda rs, pos: {
        "rsId": rs,
        "locations": [{"chromosomeName": "11", "chromosomePosition": pos}],
        "_links": {"snpAssociations": {"href": f"https://x/assoc/{rs}{{?projection}}"}},
    }

    def page1(m, u, p, b):
        return {
            "_embedded": {"singleNucleotidePolymorphisms": [snp("rs1", 100)]},
            "_links": {"next": {"href": "https://x/region?page=1"}},
        }, None, False

    def page2(m, u, p, b):
        return {
            "_embedded": {"singleNucleotidePolymorphisms": [snp("rs3", 300)]},
            "_links": {},
        }, None, False

    def assoc(m, u, p, b):
        pv = {"rs1": 1e-9, "rs3": 1e-6}[u.split("/assoc/")[1]]
        return {
            "_embedded": {"associations": [{"pvalue": pv, "efoTraits": [{"trait": "HbF, level"}]}]}
        }, None, False

    calls = []

    def mock_fetch(method, url, params=None, json_body=None):
        calls.append((method, url))
        if "findByChromBpLocationRange" in url:
            return page1(method, url, params, json_body)
        if "region?page=1" in url:
            return page2(method, url, params, json_body)
        if "/assoc/" in url:
            return assoc(method, url, params, json_body)
        return None, "http_404", False

    gwas = collect_gwas_evidence("11", [100, 300, 500], fetch=mock_fetch, window_bp=25000)
    assert set(gwas["by_pos"]) == {100, 300}
    # Templated hrefs must be cleaned
    assert not any("{" in c[1] for c in calls)

    recs = build_gwas_records(gwas, variants, p_threshold=5e-8)
    r100 = [r for r in recs if r.pos == 100][0]
    r300 = [r for r in recs if r.pos == 300][0]
    r500 = [r for r in recs if r.pos == 500][0]

    assert r100.missingness_flag == "observed_positive"
    assert r300.missingness_flag == "observed_null"
    assert json.loads(r300.evidence_value)["reason"] == "subthreshold_association"
    assert r500.missingness_flag == "absent"
    assert json.loads(r500.evidence_value)["reason"] == "not_in_catalog"
