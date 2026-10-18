"""GWAS Catalog evidence collection for Layer 1.

The catalog only stores associations that cleared a reporting threshold, so a
variant with no catalog record was never shown to have been tested -> ``absent``
(``not_in_catalog``).  A catalog record whose associations are all above
``gwas_p_threshold`` proves the variant *was* tested -> ``observed_null``
(``subthreshold_association``).

Key design decisions
--------------------
- Position-only join: the catalog doesn't report alleles in a directly
  comparable form.  A hit means the candidate *is itself* a catalogued lead
  SNP (no LD proxies).
- Windowed queries: candidates are clustered using
  :func:`~regvar.layer1.variants.plan_windows` (default 25 kb) so only
  candidate-dense sub-regions are queried, avoiding wide-locus timeouts.
- Per-window error tracking: a failure in one window marks only its enclosed
  variants as ``incomplete_region_query``; healthy windows are unaffected.
- HAL pagination: ``_links.next`` is followed; templated hrefs are cleaned.
- Association follow-up only fires for SNPs at candidate positions.
"""
from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Set

from regvar.layer1.schema import make_record
from regvar.layer1.variants import Variant, plan_windows


GWAS_BASE = "https://www.ebi.ac.uk/gwas/rest/api"


# ---------------------------------------------------------------------------
# Internal utilities
# ---------------------------------------------------------------------------


def _hal_embedded(payload: Optional[dict], key: str) -> list:
    return ((payload or {}).get("_embedded") or {}).get(key) or []


def _clean_href(href: Optional[str]) -> Optional[str]:
    """Strip HAL URL templates such as ``{?projection}``."""
    return href.split("{")[0] if href else href


def _fetch_paginated_hal(
    url: str,
    params,
    embedded_key: str,
    fetch: Callable,
    max_pages: int = 200,
):
    """Follow HAL ``_links.next`` pagination.

    Returns ``(items, error)``; ``error is not None`` means the list is
    incomplete (e.g. hit ``max_pages`` or encountered a network error).
    """
    items: list = []
    next_url: Optional[str] = url
    next_params = dict(params or {})
    pages = 0
    while next_url and pages < max_pages:
        payload, err, _ = fetch("GET", next_url, params=next_params)
        if err:
            return items, err
        items.extend(_hal_embedded(payload, embedded_key))
        nxt = _clean_href(
            ((payload or {}).get("_links") or {}).get("next", {}).get("href")
        )
        next_url, next_params = nxt, None  # next href already carries its query string
        pages += 1
    if next_url:
        return items, f"truncated_after_{max_pages}_pages"
    return items, None


def _gwas_snp_position(snp: dict, chrom: str) -> Optional[int]:
    for loc in snp.get("locations") or []:
        if (
            str(loc.get("chromosomeName")) == str(chrom)
            and loc.get("chromosomePosition") is not None
        ):
            return int(loc["chromosomePosition"])
    return None


def _assoc_pvalue(assoc: dict) -> Optional[float]:
    p = assoc.get("pvalue")
    if p is not None:
        try:
            return float(p)
        except (TypeError, ValueError):
            pass
    m, e = assoc.get("pvalueMantissa"), assoc.get("pvalueExponent")
    if m is not None and e is not None:
        try:
            return float(m) * 10 ** float(e)
        except (TypeError, ValueError):
            pass
    return None


def _resolve_traits(assoc: dict, fetch: Callable) -> List[str]:
    inline = assoc.get("efoTraits")
    if inline:
        return [t.get("trait") for t in inline if t.get("trait")]
    emb = _hal_embedded(assoc, "efoTraits")
    if emb:
        return [t.get("trait") for t in emb if t.get("trait")]
    href = _clean_href(
        (((assoc.get("_links") or {}).get("efoTraits")) or {}).get("href")
    )
    if href:
        payload, err, _ = fetch("GET", href)
        if not err:
            return [t.get("trait") for t in _hal_embedded(payload, "efoTraits") if t.get("trait")]
    return []


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def collect_gwas_evidence(
    chrom: str,
    positions: List[int],
    fetch: Callable,
    size: int = 500,
    window_bp: int = 25_000,
) -> Dict[str, Any]:
    """Query GWAS Catalog for all candidate positions on *chrom*.

    Returns a result dict::

        {
            "windows": [{"chrom", "start", "stop", "error", "n_snps"}, ...],
            "by_pos":  {pos: {"rs_ids", "assocs", "error"}, ...},
            "n_snps_in_region": int,
        }
    """
    result: Dict[str, Any] = {"windows": [], "by_pos": {}, "n_snps_in_region": 0}
    if not positions:
        return result

    wanted: Set[int] = set(positions)
    all_snps: list = []

    for s, e in plan_windows(positions, window_bp):
        # Pad by +-1 bp to guarantee inclusive matching at window boundaries.
        start, stop = max(1, s - 1), e + 1
        snps, err = _fetch_paginated_hal(
            f"{GWAS_BASE}/singleNucleotidePolymorphisms/search/findByChromBpLocationRange",
            {"chrom": str(chrom), "bpStart": start, "bpEnd": stop, "size": size},
            "singleNucleotidePolymorphisms",
            fetch,
        )
        win: Dict[str, Any] = {
            "chrom": str(chrom),
            "start": s,
            "stop": e,
            "error": err,
            "n_snps": len(snps),
        }
        result["windows"].append(win)
        all_snps.extend(snps)

    result["n_snps_in_region"] = len(all_snps)

    for snp in all_snps:
        pos = _gwas_snp_position(snp, chrom)
        if pos not in wanted:
            continue
        entry = result["by_pos"].setdefault(
            pos, {"rs_ids": [], "assocs": [], "error": None}
        )
        rs = snp.get("rsId")
        if rs and rs not in entry["rs_ids"]:
            entry["rs_ids"].append(rs)
        links = snp.get("_links") or {}
        href = _clean_href(
            ((links.get("snpAssociations") or links.get("associations")) or {}).get("href")
        )
        if not href:
            entry["error"] = entry["error"] or "no_association_link"
            continue
        assocs, aerr = _fetch_paginated_hal(href, None, "associations", fetch)
        if aerr:
            entry["error"] = aerr
        for a in assocs:
            entry["assocs"].append(
                {
                    "traits": _resolve_traits(a, fetch) or ["unresolved"],
                    "pvalue": _assoc_pvalue(a),
                    "rs_id": rs,
                }
            )
    return result


def build_gwas_records(
    gwas: Dict[str, Any],
    variants: List[Variant],
    p_threshold: float,
):
    """Convert the *gwas* result dict into :class:`~regvar.layer1.schema.EvidenceRecord` objects."""
    records = []
    for v in variants:
        entry = gwas["by_pos"].get(v.pos)
        # Find the window covering this variant for per-window error attribution.
        win_matches = [
            w
            for w in gwas.get("windows", [])
            if w["chrom"] == str(v.chrom) and w["start"] <= v.pos <= w["stop"]
        ]
        win_err = win_matches[0]["error"] if win_matches else None

        if entry is None:
            if win_err:
                records.append(
                    make_record(
                        v, "gwas", "experimental", "absent",
                        reason="incomplete_region_query", error=win_err,
                    )
                )
            else:
                records.append(
                    make_record(
                        v, "gwas", "experimental", "absent", reason="not_in_catalog"
                    )
                )
            continue

        ctx = (entry["rs_ids"] or [None])[0]
        sig = [a for a in entry["assocs"] if a["pvalue"] is not None and a["pvalue"] <= p_threshold]
        sub = [a for a in entry["assocs"] if a["pvalue"] is not None and a["pvalue"] > p_threshold]

        if sig:
            for a in sig:
                records.append(
                    make_record(
                        v, "gwas", "experimental", "observed_positive",
                        context=a.get("rs_id") or ctx,
                        traits=a["traits"],
                        pvalue=a["pvalue"],
                    )
                )
        elif sub:
            records.append(
                make_record(
                    v, "gwas", "experimental", "observed_null",
                    context=ctx,
                    reason="subthreshold_association",
                    best_pvalue=min(a["pvalue"] for a in sub),
                    n_assoc=len(sub),
                )
            )
        else:
            records.append(
                make_record(
                    v, "gwas", "experimental", "absent",
                    context=ctx,
                    reason="query_failed" if entry["error"] else "no_usable_association_records",
                    error=entry["error"],
                )
            )
    return records
