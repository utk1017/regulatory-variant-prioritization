"""GTEx Portal API v2 eQTL collection for Layer 1.

Testability rule
----------------
GTEx only tests variants with MAF >= ~1% (minor allele seen in >= 10 samples).
A *rare* variant with no eQTL hit was never tested, so its missingness flag is
``absent`` (``not_testable_rare_or_unknown_frequency``).  Only a *common*
variant (gnomAD MAF >= ``eqtl_testable_min_maf``, used as a proxy for the GTEx
cohort MAF) that has no significant hit is labelled ``observed_null``
(``common_variant_not_significant``).

Variant identity
----------------
GTEx ``variantId`` strings have the form ``chr11_5227002_T_C_b38``.  These are
parsed by :func:`parse_gtex_variant_id` and then normalised with
:func:`~regvar.layer1.variants.canon_key` for join-safe lookup.

Gene resolution
---------------
The GENCODE versioned ID (e.g. ``ENSG00000244734.11``) is resolved with
``gencodeVersion`` pinning so the versioned ID matches the chosen dataset
(gtex_v8 uses GENCODE v26).  Falls back to symbol-only lookup with a warning
if the versioned query finds nothing.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Any, Callable, Dict, List, Optional, Tuple

from regvar.layer1.schema import make_record
from regvar.layer1.variants import Variant, canon_key


GTEX_BASE = "https://gtexportal.org/api/v2"


# ---------------------------------------------------------------------------
# Gene resolution
# ---------------------------------------------------------------------------


def resolve_versioned_gencode_id(
    gene_symbol: str,
    ensembl_id: Optional[str],
    gencode_version: str,
    fetch: Callable,
) -> Tuple[Optional[str], Optional[str]]:
    """Return ``(gencode_id, note)``.

    *gencode_id* is ``None`` on failure; *note* describes how the match was
    made (or the error reason when ``None`` is returned).
    """
    attempts = (
        {"geneId": gene_symbol, "gencodeVersion": gencode_version, "genomeBuild": "GRCh38/hg38"},
        {"geneId": gene_symbol},
    )
    for params in attempts:
        payload, err, _ = fetch("GET", f"{GTEX_BASE}/reference/gene", params=params)
        if err:
            continue
        records = (payload or {}).get("data") or []
        # Prefer exact Ensembl ID match
        for rec in records:
            gid = rec.get("gencodeId") or ""
            if ensembl_id and gid.split(".")[0] == ensembl_id:
                return gid, None
        # Fall back to symbol match
        for rec in records:
            if rec.get("geneSymbol") == gene_symbol and rec.get("gencodeId"):
                return rec["gencodeId"], "matched_by_symbol_not_ensembl_id"
    return None, "gene_not_resolved"


# ---------------------------------------------------------------------------
# Variant ID parsing
# ---------------------------------------------------------------------------


def parse_gtex_variant_id(variant_id: Optional[str]) -> Optional[Tuple[str, int, str, str]]:
    """Parse a GTEx ``variantId`` string into ``(chrom, pos, ref, alt)``.

    Returns ``None`` if the string is malformed.
    """
    if not variant_id:
        return None
    parts = variant_id.split("_")
    if len(parts) < 4:
        return None
    try:
        pos = int(parts[1])
    except ValueError:
        return None
    return parts[0].replace("chr", ""), pos, parts[2], parts[3]


# ---------------------------------------------------------------------------
# Paginated eQTL fetch
# ---------------------------------------------------------------------------


def fetch_tissue_eqtls(
    gencode_id: str,
    tissue: str,
    dataset: str,
    fetch: Callable,
    items_per_page: int,
    max_pages: int,
):
    """Return all significant eQTL rows for one gene x tissue.

    Returns ``(rows, error)``; *error* is ``None`` on a clean finish.
    """
    rows: list = []
    seen: set = set()
    for page in range(max_pages):
        payload, err, _ = fetch(
            "GET",
            f"{GTEX_BASE}/association/singleTissueEqtl",
            params={
                "gencodeId": gencode_id,
                "tissueSiteDetailId": tissue,
                "datasetId": dataset,
                "itemsPerPage": items_per_page,
                "page": page,
            },
        )
        if err:
            return rows, err
        data = (payload or {}).get("data") or []
        new = 0
        for rec in data:
            vid = rec.get("variantId")
            if vid in seen:
                continue
            seen.add(vid)
            rows.append(rec)
            new += 1
        n_pages = ((payload or {}).get("paging_info") or {}).get("numberOfPages")
        if len(data) < items_per_page:
            return rows, None
        if new == 0:
            return rows, "pagination_not_advancing"
        if n_pages is not None and page + 1 >= n_pages:
            return rows, None
    return rows, f"truncated_after_{max_pages}_pages"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def collect_eqtl_evidence(
    gene_meta: Dict[str, Any],
    tissues: List[str],
    dataset: str,
    gencode_version: str,
    items_per_page: int,
    max_pages: int,
    fetch: Callable,
) -> Dict[str, Any]:
    """Query GTEx for all configured tissues.

    Returns::

        {
            "gencode_id": str or None,
            "resolve_note": str or None,
            "resolve_error": str or None,
            "tissues": {
                tissue: {"error": ..., "n_rows": ..., "by_key": {canon_key: [hits]}},
            },
        }
    """
    out: Dict[str, Any] = {
        "gencode_id": None,
        "resolve_note": None,
        "resolve_error": None,
        "tissues": {},
    }
    gid, note = resolve_versioned_gencode_id(
        gene_meta["gene_symbol"], gene_meta.get("gene_id"), gencode_version, fetch
    )
    if gid is None:
        out["resolve_error"] = note
        return out
    out["gencode_id"], out["resolve_note"] = gid, note

    for tissue in tissues:
        rows, err = fetch_tissue_eqtls(gid, tissue, dataset, fetch, items_per_page, max_pages)
        by_key: dict = defaultdict(list)
        for rec in rows:
            parsed = parse_gtex_variant_id(rec.get("variantId"))
            if not parsed:
                continue
            by_key[canon_key(*parsed)].append(
                {
                    "slope": rec.get("nes", rec.get("slope")),
                    "pvalue": rec.get("pValue", rec.get("pval")),
                    "variant_id": rec.get("variantId"),
                }
            )
        out["tissues"][tissue] = {"error": err, "n_rows": len(rows), "by_key": dict(by_key)}
    return out


def build_eqtl_records(
    eqtl: Dict[str, Any],
    variants: List[Variant],
    tissues: List[str],
    maf_by_key: Dict,
    min_maf: float,
):
    """Convert the *eqtl* result dict into :class:`~regvar.layer1.schema.EvidenceRecord` objects."""
    records = []
    for tissue in tissues:
        tinfo = eqtl["tissues"].get(tissue)
        for v in variants:
            if eqtl.get("resolve_error") or tinfo is None:
                records.append(
                    make_record(
                        v, "eqtl", "experimental", "absent",
                        context=tissue,
                        reason="query_failed",
                        error=eqtl.get("resolve_error") or "tissue_not_queried",
                    )
                )
                continue
            hits = tinfo["by_key"].get(v.key)
            if hits:
                for h in hits:
                    records.append(
                        make_record(
                            v, "eqtl", "experimental", "observed_positive",
                            context=tissue,
                            slope=h["slope"],
                            pvalue=h["pvalue"],
                            gtex_variant_id=h["variant_id"],
                        )
                    )
                continue
            if tinfo["error"]:
                records.append(
                    make_record(
                        v, "eqtl", "experimental", "absent",
                        context=tissue,
                        reason="incomplete_tissue_query",
                        error=tinfo["error"],
                    )
                )
                continue
            maf = maf_by_key.get(v.key)
            if maf is not None and maf >= min_maf:
                records.append(
                    make_record(
                        v, "eqtl", "experimental", "observed_null",
                        context=tissue,
                        reason="common_variant_not_significant",
                        maf=maf,
                    )
                )
            else:
                records.append(
                    make_record(
                        v, "eqtl", "experimental", "absent",
                        context=tissue,
                        reason="not_testable_rare_or_unknown_frequency",
                        maf=maf,
                    )
                )
    return records
