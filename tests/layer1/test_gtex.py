import json
import pandas as pd
from regvar.layer1.gtex import (
    build_eqtl_records,
    collect_eqtl_evidence,
    fetch_tissue_eqtls,
    resolve_versioned_gencode_id,
)
from regvar.layer1.variants import Variant, canon_key, to_variants


def _variants(rows):
    return to_variants(pd.DataFrame([{"chrom": "11", "pos": p, "ref": r, "alt": a} for p, r, a in rows]))


def test_resolve_versioned_gencode_id():
    def mock_fetch(method, url, params=None, json_body=None):
        if params and "gencodeVersion" in params:
            return {
                "data": [
                    {"geneSymbol": "HBB", "gencodeId": "ENSG00000244734.11"},
                    {"geneSymbol": "OTHER", "gencodeId": "ENSG00000999999.1"},
                ]
            }, None, False
        return None, "http_404", False

    gid, note = resolve_versioned_gencode_id("HBB", "ENSG00000244734", "v26", mock_fetch)
    assert gid == "ENSG00000244734.11"
    assert note is None


def test_gtex_testability_and_failed_tissue():
    variants = _variants([(100, "A", "G"), (200, "C", "T"), (300, "G", "A"), (400, "T", "C")])
    eqtl = {
        "gencode_id": "ENSG1.1",
        "resolve_note": None,
        "resolve_error": None,
        "tissues": {
            "Whole_Blood": {
                "error": None,
                "n_rows": 1,
                "by_key": {
                    canon_key("11", 100, "A", "G"): [
                        {"slope": 0.4, "pvalue": 1e-8, "variant_id": "chr11_100_A_G_b38"}
                    ]
                },
            },
            "Liver": {"error": "http_500", "n_rows": 0, "by_key": {}},
        },
    }
    maf = {canon_key("11", 200, "C", "T"): 0.2, canon_key("11", 300, "G", "A"): 0.0001}

    recs = build_eqtl_records(eqtl, variants, ["Whole_Blood", "Liver"], maf, min_maf=0.01)
    by_pos_tissue = {(r.pos, r.context): r for r in recs}

    # Pos 100: Hit in Whole_Blood
    assert by_pos_tissue[(100, "Whole_Blood")].missingness_flag == "observed_positive"

    # Pos 200: Common non-hit (MAF 0.2 >= 0.01) -> absent / unverified (not tested-null)
    assert by_pos_tissue[(200, "Whole_Blood")].missingness_flag == "absent"
    pos200_val = json.loads(by_pos_tissue[(200, "Whole_Blood")].evidence_value)
    assert pos200_val["reason"] == "no_significant_eqtl_association_returned_test_status_unverified"
    assert pos200_val["test_status"] == "unverified"

    # Pos 300: Rare variant (MAF 0.0001 < 0.01) -> absent (not testable)
    assert by_pos_tissue[(300, "Whole_Blood")].missingness_flag == "absent"
    assert json.loads(by_pos_tissue[(300, "Whole_Blood")].evidence_value)["reason"].startswith("not_testable")

    # Pos 400: Not in gnomAD at all -> absent
    assert by_pos_tissue[(400, "Whole_Blood")].missingness_flag == "absent"

    # All Liver queries failed -> absent
    assert all(r.missingness_flag == "absent" for r in recs if r.context == "Liver")


def test_empty_gtex_result_for_common_variant_is_not_a_tested_null():
    """An empty singleTissueEqtl query cannot establish a tested-negative result."""
    v = Variant("11", 900, "A", "G", canon_key("11", 900, "A", "G"))
    for gene_rows in (0, 3):
        eqtl_fixture = {
            "gencode_id": "ENSG00000244734.3",
            "resolve_error": None,
            "tissues": {
                "Whole_Blood": {"error": None, "n_rows": gene_rows, "by_key": {}},
            },
        }
        records = build_eqtl_records(
            eqtl_fixture, [v], ["Whole_Blood"], {v.key: 0.12}, 0.01
        )
        assert len(records) == 1
        assert records[0].missingness_flag == "absent"
        payload = json.loads(records[0].evidence_value)
        assert payload["reason"] == "no_significant_eqtl_association_returned_test_status_unverified"
        assert payload["test_status"] == "unverified"


def test_fetch_tissue_eqtls_pagination():
    rows = [{"variantId": f"chr11_{i}_A_G_b38", "nes": 0.1, "pValue": 1e-6} for i in range(5)]

    def mock_paged(method, url, params=None, json_body=None):
        page = params["page"]
        chunk = rows[page * 2 : (page + 1) * 2]
        return {"data": chunk, "paging_info": {"numberOfPages": 3, "page": page}}, None, False

    got, err = fetch_tissue_eqtls("E.1", "Whole_Blood", "gtex_v8", mock_paged, items_per_page=2, max_pages=10)
    assert err is None
    assert len(got) == 5
