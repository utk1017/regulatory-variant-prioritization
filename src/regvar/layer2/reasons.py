"""Reason classification and multi-context aggregation rules for Layer 2.

Handles parsing evidence payloads, classifying absent reasons (biological gap vs.
experimental null vs. transient error vs. unknown reason), aggregating status
across contexts (positive > null > absent), and scoping records to target tissues/cells.
"""
from __future__ import annotations

import json
from typing import Any, Dict, List, Optional, Sequence, Tuple

from regvar.layer2.schema import REASON_CATEGORIES


def _norm_ctx(ctx: Any) -> str:
    """Normalize context values: None, NaN, empty strings map to 'default'."""
    if ctx is None or (isinstance(ctx, float) and ctx != ctx):
        return "default"
    s = str(ctx).strip()
    return "default" if not s or s.lower() == "nan" else s


def parse_evidence_value(val: Any) -> Dict[str, Any]:
    """Safely parse a JSON evidence_value payload or return empty dict."""
    if isinstance(val, dict):
        return val
    if not isinstance(val, str) or not val.strip():
        return {}
    try:
        parsed = json.loads(val)
        return parsed if isinstance(parsed, dict) else {"raw": parsed}
    except Exception:
        return {"raw": str(val)}


def classify_evidence_reason(
    missingness_flag: str,
    reason: Optional[str] = None,
    error: Optional[str] = None,
) -> str:
    """Classify an evidence observation into one of:
    - POSITIVE (observed_positive hit)
    - NULL_EXPERIMENTAL (observed negative result)
    - GAP_BIOLOGICAL (true missingness eligible for computational fill)
    - ERROR_TRANSIENT (network/query failure requiring pipeline rerun)
    - UNKNOWN_REASON (unrecognized reason string treated conservatively)
    """
    if missingness_flag == "observed_positive":
        return "POSITIVE"
    if missingness_flag == "observed_null":
        return "NULL_EXPERIMENTAL"

    # For absent observations:
    if error or reason in ("query_failed", "incomplete_region_query", "incomplete_tissue_query"):
        return "ERROR_TRANSIENT"
    if reason is None:
        return "GAP_BIOLOGICAL"
    if reason in REASON_CATEGORIES:
        return REASON_CATEGORIES[reason]
    return "UNKNOWN_REASON"


def aggregate_modality_status(
    records: Sequence[Dict[str, Any]]
) -> Tuple[str, str, List[str], Dict[str, Dict[str, Any]]]:
    """Aggregate multiple EvidenceRecord rows for one (variant, modality) across contexts.

    Precedence:
    1. 'observed_positive' if ANY record is observed_positive
    2. 'observed_null' if no positive, but AT LEAST ONE record is observed_null
    3. 'absent' if ALL records are absent

    Returns:
        (aggregate_status, primary_category, list_of_reasons, dict_by_context)
    """
    if not records:
        return "absent", "GAP_BIOLOGICAL", ["untested"], {}

    by_context: Dict[str, Dict[str, Any]] = {}
    flags = [r["missingness_flag"] for r in records]
    reasons: List[str] = []
    errors: List[str] = []

    for r in records:
        ctx = _norm_ctx(r.get("context"))
        ev = parse_evidence_value(r.get("evidence_value"))
        reason = ev.get("reason")
        err = ev.get("error")
        cat = classify_evidence_reason(r["missingness_flag"], reason, err)
        by_context[ctx] = {
            "flag": r["missingness_flag"],
            "reason": reason,
            "category": cat,
            "error": err,
            "evidence_value": ev,
        }
        if reason:
            reasons.append(reason)
        if err:
            errors.append(err)

    if "observed_positive" in flags:
        status = "observed_positive"
        category = "POSITIVE"
    elif "observed_null" in flags:
        status = "observed_null"
        category = "NULL_EXPERIMENTAL"
    else:
        status = "absent"
        context_cats = {c["category"] for c in by_context.values()}
        if errors or "ERROR_TRANSIENT" in context_cats:
            category = "ERROR_TRANSIENT"
        elif "UNKNOWN_REASON" in context_cats:
            category = "UNKNOWN_REASON"
        else:
            category = "GAP_BIOLOGICAL"

    return status, category, reasons, by_context


def filter_records_to_targets(
    records: Sequence[Dict[str, Any]],
    target_tissues: Optional[Sequence[str]] = None,
    target_cell_types: Optional[Sequence[str]] = None,
) -> List[Dict[str, Any]]:
    """Filter records to target tissues (for eQTL) and target cell types (for MPRA/CRISPR).

    Default/context-less records (None, NaN, 'default') are preserved as they represent
    global element or library properties.
    """
    kept: List[Dict[str, Any]] = []
    for r in records:
        ctx = _norm_ctx(r.get("context"))
        mod = r["modality"]
        if ctx != "default":
            if mod == "eqtl" and target_tissues is not None and ctx not in target_tissues:
                continue
            if mod in ("mpra", "crispr") and target_cell_types is not None and ctx not in target_cell_types:
                continue
        kept.append(r)
    return kept
