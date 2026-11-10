"""GWAS Catalog evidence collection for Layer 1.

The catalog only stores associations that cleared a reporting threshold, so a
variant with no catalog record was never shown to have been tested -> absent
(not_in_catalog). A catalog record whose associations are all above
gwas_p_threshold proves the variant *was* tested -> observed_null
(subthreshold_association).

Supports:
- Auto fallback between REST API and official bulk release TSV/zip (ftp.ebi.ac.uk)
- Diagnostic endpoint probing for 404/410 errors
- Explicit position-level vs. allele-specific matching markers
"""
from __future__ import annotations

import os
import re
import time
import zipfile
from typing import Any, Callable, Dict, List, Optional, Set

import pandas as pd
import requests

from regvar.common.provenance import _file_hash
from regvar.layer1.schema import make_record
from regvar.layer1.variants import Variant, plan_windows

GWAS_BASE = "https://www.ebi.ac.uk/gwas/rest/api/v2"


def _hal_embedded(payload: Optional[dict | list], key: str) -> list:
    """Item list from a HAL (_embedded[key]) or plain (key / content / results / data / items) payload."""
    if isinstance(payload, list):
        return payload
    if not isinstance(payload, dict):
        return []
    emb = payload.get("_embedded")
    if isinstance(emb, dict) and isinstance(emb.get(key), list):
        return emb[key]
    for k in (key, "content", "results", "data", "items"):
        if isinstance(payload.get(k), list):
            return payload[k]
    return []


def _unrecognized_payload(payload: Any, key: str) -> bool:
    """True when a 200 response is a dict with neither the expected item container nor any HAL/paging
    marker -- the API shape changed. That must be an ERROR, never silently read as 'no results'."""
    if payload is None or isinstance(payload, list):
        return False
    if not isinstance(payload, dict):
        return True
    markers = ("_embedded", "_links", "links", "page", key, "content", "results", "data", "items")
    return bool(payload) and not any(m in payload for m in markers)


def _next_href(payload: Any) -> Optional[str]:
    if not isinstance(payload, dict):
        return None
    for container in ("_links", "links"):
        links = payload.get(container)
        nxt = links.get("next") if isinstance(links, dict) else None
        if isinstance(nxt, dict) and nxt.get("href"):
            return nxt["href"]
        if isinstance(nxt, str) and nxt:
            return nxt
    nxt = payload.get("next")
    return nxt if isinstance(nxt, str) and nxt else None


def _clean_href(href: Optional[str]) -> Optional[str]:
    return href.split("{")[0] if href else href


def fetch_paginated_hal(
    url: str,
    params: Optional[dict],
    embedded_key: str,
    fetch: Callable,
    max_pages: int = 200,
):
    """Follow HAL _links.next and collect all pages."""
    items = []
    next_url = url
    next_params = dict(params or {})
    pages = 0

    while next_url and pages < max_pages:
        payload, err, _ = fetch("GET", next_url, params=next_params)
        if err:
            return items, err

        if _unrecognized_payload(payload, embedded_key):
            return items, "unrecognized_response_shape"
        items.extend(_hal_embedded(payload, embedded_key))

        next_url = _clean_href(_next_href(payload))
        next_params = None
        pages += 1

    if next_url:
        return items, f"truncated_after_{max_pages}_pages"

    return items, None


def gwas_snp_position(snp: dict, chrom: str) -> Optional[int]:
    for loc in snp.get("locations") or []:
        if (
            str(loc.get("chromosome_name")) == str(chrom)
            and loc.get("chromosome_position") is not None
        ):
            return int(loc["chromosome_position"])
    return None


def assoc_pvalue(assoc: dict) -> Optional[float]:
    p = assoc.get("p_value")
    if p is not None:
        try:
            return float(p)
        except (TypeError, ValueError):
            pass

    m = assoc.get("pvalue_mantissa")
    e = assoc.get("pvalue_exponent")

    if m is not None and e is not None:
        try:
            return float(m) * 10 ** float(e)
        except (TypeError, ValueError):
            pass

    return None


def resolve_traits(assoc: dict, fetch: Optional[Callable] = None) -> List[str]:
    traits = assoc.get("efo_traits") or []
    names = [
        t.get("efo_trait")
        for t in traits
        if isinstance(t, dict) and t.get("efo_trait")
    ]

    if not names:
        names = [
            t for t in (assoc.get("reported_trait") or [])
            if isinstance(t, str) and t
        ]

    return names


def diagnose_gwas_endpoint(chrom: str, start: int, stop: int, timeout: int = 20):
    """One raw probe of the SNP-by-region call. Runs automatically when a GWAS window fails and prints
    what the server actually returned (status, final URL, content type, top-level JSON keys)."""
    url = f"{GWAS_BASE}/single-nucleotide-polymorphisms"
    params = {"chromosome": str(chrom), "bp_start": int(start), "bp_end": int(stop), "size": 5, "page": 0}
    try:
        r = requests.get(url, params=params, timeout=timeout, headers={"Accept": "application/json"})
    except requests.exceptions.RequestException as e:
        print(f"GWAS diagnostic probe: request failed ({type(e).__name__}: {e})")
        return None
    try:
        j = r.json()
        keys = sorted(j.keys()) if isinstance(j, dict) else f"list[{len(j)}]"
    except ValueError:
        keys = None
    info = {
        "status": r.status_code,
        "url": r.url,
        "content_type": r.headers.get("Content-Type"),
        "top_level_keys": keys,
        "body_head": r.text[:200],
    }
    print("GWAS diagnostic probe:")
    for k, v in info.items():
        print(f"  {k}: {v}")
    if r.status_code in (404, 410):
        print(
            "  -> The endpoint/path is gone or renamed. Open https://www.ebi.ac.uk/gwas/rest/api/v2/swagger-ui/index.html, "
            "find the SNP/variant region filter, and/or keep layer1.gwas_source: auto so the official bulk release file is used instead."
        )
    return info


def collect_gwas_evidence(chrom: str, positions: Sequence[int], fetch: Callable, size: int = 500, window_bp: int = 25000):
    result = {"windows": [], "by_pos": {}, "n_snps_in_region": 0}

    if not positions:
        return result

    wanted = set(positions)
    windows = plan_windows(positions, window_bp)
    all_snps = []

    for s, e in windows:
        start, stop = max(1, s - 1), e + 1

        snps, err = fetch_paginated_hal(
            f"{GWAS_BASE}/single-nucleotide-polymorphisms",
            {
                "chromosome": str(chrom),
                "bp_start": start,
                "bp_end": stop,
                "size": size,
                "page": 0,
            },
            "snps",
            fetch,
        )

        result["windows"].append({
            "chrom": str(chrom),
            "start": s,
            "stop": e,
            "error": err,
            "n_snps": len(snps),
        })

        all_snps.extend(snps)

    # Deduplicate SNPs because adjacent windows can overlap.
    unique_snps = {}
    for snp in all_snps:
        rs = snp.get("rs_id")
        pos = gwas_snp_position(snp, chrom)
        key = (rs, pos) if rs else (None, pos, str(snp))
        unique_snps[key] = snp

    all_snps = list(unique_snps.values())
    result["n_snps_in_region"] = len(all_snps)

    for snp in all_snps:
        pos = gwas_snp_position(snp, chrom)

        if pos not in wanted:
            continue

        entry = result["by_pos"].setdefault(
            pos, {"rs_ids": [], "assocs": [], "error": None}
        )

        rs = snp.get("rs_id")
        if rs and rs not in entry["rs_ids"]:
            entry["rs_ids"].append(rs)

        if not rs:
            entry["error"] = entry["error"] or "missing_rs_id"
            continue

        assocs, aerr = fetch_paginated_hal(
            f"{GWAS_BASE}/associations",
            {"rs_id": rs, "size": size, "page": 0},
            "associations",
            fetch,
        )

        if aerr:
            entry["error"] = aerr

        for assoc in assocs:
            entry["assocs"].append({
                "traits": resolve_traits(assoc) or ["unresolved"],
                "pvalue": assoc_pvalue(assoc),
                "rs_id": rs,
            })

    return result


def _gwas_norm_col(c: Any) -> str:
    return re.sub(r"\s+", "_", str(c).strip().upper())


_GWAS_WANTED = {"CHR_ID", "CHR_POS", "SNPS", "P-VALUE", "PVALUE_MLOG", "DISEASE/TRAIT", "MAPPED_TRAIT"}


def _split_multi(value: Any) -> List[str]:
    """CHR_ID / CHR_POS / SNPS hold several values for multi-SNP rows ('a; b' or 'a x b')."""
    if not isinstance(value, str):
        return []
    return [t.strip() for t in re.split(r"\s*;\s*|\s+x\s+", value) if t.strip()]


def download_gwas_catalog(url: str, cache_dir: str, ttl_days: int = 14, timeout: int = 300) -> Tuple[Optional[str], Optional[str]]:
    """Download (or reuse) the official GWAS Catalog bulk associations file. Returns (path, error).
    A stale cached copy is reused if the download fails (its age is recorded by the caller)."""
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, os.path.basename(url))
    fresh = os.path.exists(path) and (time.time() - os.path.getmtime(path)) < ttl_days * 86400
    if fresh:
        return path, None
    tmp = path + ".part"
    try:
        with requests.get(url, stream=True, timeout=timeout) as r:
            r.raise_for_status()
            with open(tmp, "wb") as fh:
                for chunk in r.iter_content(1 << 20):
                    fh.write(chunk)
        os.replace(tmp, path)
        return path, None
    except (requests.exceptions.RequestException, OSError) as e:
        if os.path.exists(path):
            print(f"WARNING: bulk GWAS download failed ({type(e).__name__}); using older cached copy.")
            return path, None
        return None, f"bulk_download_failed:{type(e).__name__}"


def collect_gwas_evidence_from_download(chrom: str | int, positions: Sequence[int], path: str, window_bp: int = 25000):
    """Same output shape as collect_gwas_evidence(), read from the official bulk TSV (zip).
    Still position-level (the file lists the strongest SNP per association), not allele-specific."""
    result = {
        "windows": [],
        "by_pos": {},
        "n_snps_in_region": 0,
        "source_label": "GWAS Catalog bulk download (associations file)",
    }
    if not positions:
        return result
    wanted = set(positions)
    windows = plan_windows(positions, window_bp)
    chrom = str(chrom)
    frames = []
    try:
        with zipfile.ZipFile(path) as zf:
            members = [m for m in zf.namelist() if m.lower().endswith((".tsv", ".txt"))]
            if not members:
                raise ValueError("no tsv member in zip")
            with zf.open(members[0]) as fh:
                for chunk in pd.read_csv(
                    fh,
                    sep="\t",
                    dtype=str,
                    chunksize=200_000,
                    low_memory=False,
                    usecols=lambda c: _gwas_norm_col(c) in _GWAS_WANTED,
                    on_bad_lines="skip",
                ):
                    chunk.columns = [_gwas_norm_col(c) for c in chunk.columns]
                    keep = chunk["CHR_ID"].fillna("").str.contains(rf"(?:^|\s|;){re.escape(chrom)}(?:$|\s|;)", regex=True)
                    if keep.any():
                        frames.append(chunk[keep])
    except (OSError, ValueError, KeyError, zipfile.BadZipFile, pd.errors.ParserError) as e:
        err = f"bulk_file_unreadable:{type(e).__name__}"
        result["windows"] = [{"chrom": chrom, "start": s, "stop": e_, "error": err, "n_snps": 0} for s, e_ in windows]
        return result

    snps = {}   # (rs_id, pos) -> list of association dicts
    for df in frames:
        for rec in df.to_dict("records"):
            chrs, poss, rss = (_split_multi(rec.get(k)) for k in ("CHR_ID", "CHR_POS", "SNPS"))
            n = max(len(chrs), len(poss), len(rss), 1)
            chrs = chrs * n if len(chrs) == 1 else chrs
            rss = rss + [None] * (n - len(rss))
            try:
                p = float(rec.get("P-VALUE"))
            except (TypeError, ValueError):
                p = None
            if p is None or p != p:          # missing / NaN -> derive from -log10(p) when available
                try:
                    p = 10 ** (-float(rec.get("PVALUE_MLOG")))
                except (TypeError, ValueError):
                    p = None
                if p is not None and p != p:
                    p = None
            trait = rec.get("MAPPED_TRAIT") or rec.get("DISEASE/TRAIT")
            trait = trait if isinstance(trait, str) and trait else "unresolved"
            for c, ps, rs in zip(chrs, poss, rss):
                if c != chrom or not str(ps).isdigit():
                    continue
                snps.setdefault((rs, int(ps)), []).append({"traits": [trait], "pvalue": p, "rs_id": rs})

    lo, hi = min(s for s, _ in windows) - 1, max(e for _, e in windows) + 1
    result["n_snps_in_region"] = len({k for k in snps if lo <= k[1] <= hi})
    for s, e_ in windows:
        n_here = len({k for k in snps if s - 1 <= k[1] <= e_ + 1})
        result["windows"].append({"chrom": chrom, "start": s, "stop": e_, "error": None, "n_snps": n_here})
    for (rs, pos), assocs in snps.items():
        if pos in wanted:
            entry = result["by_pos"].setdefault(pos, {"rs_ids": [], "assocs": [], "error": None})
            if rs and rs not in entry["rs_ids"]:
                entry["rs_ids"].append(rs)
            entry["assocs"].extend(assocs)
    return result


def gwas_via_bulk_file(chrom: str | int, positions: Sequence[int], L1: Dict[str, Any]):
    """Download/reuse the bulk file and extract GWAS evidence. Returns (result, file_info, error)."""
    import datetime as _dt
    cache_dir = os.path.join(os.path.dirname(str(L1.get("http_cache_dir", "data/cache/http")).rstrip("/")) or ".", "gwas_bulk")
    url = L1.get("gwas_download_url", "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip")
    ttl = L1.get("gwas_download_ttl_days", 14)
    path, err = download_gwas_catalog(url, cache_dir, ttl)
    if err:
        return None, None, err
    res = collect_gwas_evidence_from_download(chrom, positions, path, window_bp=L1.get("gwas_window_bp", 25000))
    info = {
        "url": url,
        "path": path,
        "sha256": _file_hash(str(path)),
        "file_mtime_utc": _dt.datetime.fromtimestamp(os.path.getmtime(path), _dt.timezone.utc).isoformat(),
    }
    return res, info, None


def _gwas_problems(g: Dict[str, Any]):
    errs = sorted({str(w["error"]) for w in g.get("windows", []) if w.get("error")})
    pos_errs = [p for p, e in g.get("by_pos", {}).items() if e.get("error")]
    return errs, pos_errs


def _gwas_warnings(g: Dict[str, Any], label: str):
    errs, pos_errs = _gwas_problems(g)
    out = []
    if errs:
        bad = sum(1 for w in g.get("windows", []) if w.get("error"))
        out.append(f"GWAS ({label}): {bad} of {len(g.get('windows', []))} window(s) failed query ({', '.join(errs)})")
    if pos_errs:
        out.append(f"GWAS ({label}): association query errors at {len(pos_errs)} candidate position(s)")
    return out


def collect_gwas_with_fallback(chrom: str | int, positions: Sequence[int], fetch: Callable, L1: Dict[str, Any], bulk_fn=None, probe_fn=None):
    """layer1.gwas_source: 'api' | 'download' | 'auto' (REST API first, official bulk file if it fails).
    Returns {'gwas', 'used', 'file_info', 'warnings', 'notes'}; only unresolved failures become warnings."""
    bulk_fn = bulk_fn or gwas_via_bulk_file
    probe_fn = probe_fn or diagnose_gwas_endpoint
    mode = str(L1.get("gwas_source", "auto")).lower()
    if mode not in ("auto", "api", "download"):
        raise ValueError(f"layer1.gwas_source must be auto, api or download (got {mode!r})")
    out = {"gwas": None, "used": None, "file_info": None, "warnings": [], "notes": []}
    api_bad, api_errs = False, []
    if mode in ("auto", "api"):
        api = collect_gwas_evidence(chrom, positions, fetch=fetch, window_bp=L1.get("gwas_window_bp", 25000))
        out["gwas"], out["used"] = api, "api"
        api_errs, api_pos_errs = _gwas_problems(api)
        api_bad = bool(api_errs or api_pos_errs)
        if api_errs:
            w0 = next(w for w in api["windows"] if w["error"])
            probe_fn(w0["chrom"], w0["start"], w0["stop"])
        if mode == "api" or not api_bad:
            out["warnings"] = _gwas_warnings(api, "REST API")
            return out
    bulk, info, derr = bulk_fn(chrom, positions, L1)
    bulk_errs = _gwas_problems(bulk)[0] if bulk else []
    if not derr and not bulk_errs:
        out.update(gwas=bulk, used="download", file_info=info)
        if api_bad:
            out["notes"].append(
                f"GWAS REST API failed ({', '.join(api_errs) or 'position-level errors'}); evidence taken from "
                f"the official bulk release file instead (recorded in the manifest)."
            )
        return out
    detail = derr or ", ".join(bulk_errs)
    if out["gwas"] is None:   # download-only mode, and it failed
        out["gwas"] = bulk if bulk else {
            "windows": [{"chrom": str(chrom), "start": min(positions) if positions else 0, "stop": max(positions) if positions else 0,
                         "error": derr, "n_snps": 0}],
            "by_pos": {},
            "n_snps_in_region": 0,
            "source_label": "GWAS Catalog bulk download (associations file)",
        }
        out["used"] = "download"
        out["warnings"] = _gwas_warnings(out["gwas"], "bulk file")
    else:
        out["warnings"] = _gwas_warnings(out["gwas"], "REST API") + [f"GWAS: bulk release file also unavailable ({detail})"]
    return out


def build_gwas_records(gwas: Dict[str, Any], variants: List[Variant], p_threshold: float):
    """Create one evidence record per candidate variant from GWAS Catalog position matches.

    The API query identifies SNPs by chromosome/position, not by the candidate's ref/alt alleles.
    Therefore all emitted records explicitly mark the match as position-level, not allele-specific.
    """
    records = []
    src_label = gwas.get("source_label") or "GWAS Catalog REST API v2"
    by_pos = gwas.get("by_pos") or {}
    windows = gwas.get("windows") or []
    for v in variants:
        entry = by_pos.get(v.pos)
        if entry:
            assocs = entry.get("assocs") or []
            details = []
            valid_pvalues = []
            for assoc in assocs:
                p = assoc.get("pvalue")
                if p is not None:
                    try:
                        p = float(p)
                        valid_pvalues.append(p)
                    except (TypeError, ValueError):
                        p = None
                details.append({
                    "rs_id": assoc.get("rs_id"),
                    "traits": assoc.get("traits") or ["unresolved"],
                    "pvalue": p,
                })
            common = {
                "matching_level": "position",
                "allele_specific": False,
                "candidate_alleles": {"ref": v.ref, "alt": v.alt},
                "candidate_position": int(v.pos),
                "rs_ids": entry.get("rs_ids") or [],
                "associations": details,
                "source": src_label,
                "p_value_threshold": float(p_threshold),
                "dataset_id": "gwas_catalog",
                "dataset_version": "v2",
            }
            if entry.get("error") and not assocs:
                records.append(
                    make_record(
                        v, "gwas", "experimental", "absent",
                        context="GWAS Catalog", reason="query_failed", error=entry["error"], **common
                    )
                )
            elif not assocs:
                records.append(
                    make_record(
                        v, "gwas", "experimental", "absent",
                        context="GWAS Catalog", reason="no_usable_association_records", **common
                    )
                )
            elif not valid_pvalues:
                records.append(
                    make_record(
                        v, "gwas", "experimental", "absent",
                        context="GWAS Catalog", reason="p_value_missing", **common
                    )
                )
            elif any(p <= p_threshold for p in valid_pvalues):
                sig_traits = sorted({
                    trait for a in details if a["pvalue"] is not None and a["pvalue"] <= p_threshold for trait in a["traits"]
                })
                records.append(
                    make_record(
                        v, "gwas", "experimental", "observed_positive",
                        context="GWAS Catalog", traits=sig_traits, min_pvalue=min(valid_pvalues), **common
                    )
                )
            else:
                records.append(
                    make_record(
                        v, "gwas", "experimental", "observed_null",
                        context="GWAS Catalog", reason="subthreshold_association",
                        traits=sorted({t for a in details for t in a["traits"]}),
                        min_pvalue=min(valid_pvalues), **common
                    )
                )
            continue

        failed = [
            w for w in windows
            if w.get("chrom") == v.chrom and int(w.get("start", 0)) <= v.pos <= int(w.get("stop", -1)) and w.get("error")
        ]
        if failed:
            records.append(
                make_record(
                    v, "gwas", "experimental", "absent",
                    context="GWAS Catalog", reason="incomplete_region_query",
                    errors=[w.get("error") for w in failed],
                    matching_level="position", allele_specific=False,
                    candidate_alleles={"ref": v.ref, "alt": v.alt},
                    source=src_label,
                    dataset_id="gwas_catalog",
                    dataset_version="v2",
                )
            )
        else:
            records.append(
                make_record(
                    v, "gwas", "experimental", "absent",
                    context="GWAS Catalog", reason="not_in_catalog",
                    matching_level="position", allele_specific=False,
                    candidate_alleles={"ref": v.ref, "alt": v.alt},
                    source=src_label,
                    dataset_id="gwas_catalog",
                    dataset_version="v2",
                )
            )
    return records
