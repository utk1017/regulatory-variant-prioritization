"""AlphaGenome Trigger and Agreement Queue Compilation.

Scans candidate variants for missing completable modalities (eQTL, MPRA, CRISPR)
and compiles AlphaGenome computational inference jobs. Also isolates observed
experimental hits into an independent agreement benchmark queue.
"""
from __future__ import annotations

from collections import defaultdict
import json
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

from regvar.layer2.reasons import (
    _norm_ctx,
    aggregate_modality_status,
    filter_records_to_targets,
    parse_evidence_value,
)
from regvar.layer2.schema import (
    MODALITY_META,
    NEVER_COMPLETE_CATEGORIES,
)


def _trigger_row(
    counter: int,
    v: Any,
    mod: str,
    ctx: str,
    gene_meta: Dict[str, Any],
    reason: Optional[str],
    source_ev: Any,
    head_mapping: Optional[Dict[str, Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    meta = (head_mapping or MODALITY_META).get(mod, {})
    return {
        "trigger_id": f"AG_TRIG_{counter:06d}",
        "chrom": v.chrom,
        "pos": v.pos,
        "ref": v.ref,
        "alt": v.alt,
        "variant_key": f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}",
        "modality": mod,
        "context": ctx,
        "target_gene": gene_meta.get("gene_symbol", "UNKNOWN"),
        "alphagenome_output_types": ";".join(meta.get("alphagenome_output_types", [])),
        "alphagenome_proxy_quality": meta.get("proxy_quality", "unknown"),
        "alphagenome_binding_status": "unverified",
        "receptive_field_bp": meta.get("receptive_field_bp", 1048576),
        "trigger_reason": reason,
        "source_evidence_value": json.dumps(
            source_ev if isinstance(source_ev, dict) else parse_evidence_value(source_ev)
        ),
    }


def _queried_contexts(
    evidence_by_key: Dict[Tuple[str, int, str, str], List[Dict[str, Any]]],
    modality: str,
) -> Set[str]:
    """Return all non-default contexts actually queried for this modality across variants."""
    return {
        _norm_ctx(r.get("context"))
        for recs in evidence_by_key.values()
        for r in recs
        if r["modality"] == modality
    } - {"default"}


def compile_alphagenome_triggers(
    variants: Sequence[Any],
    evidence_by_key: Dict[Tuple[str, int, str, str], List[Dict[str, Any]]],
    completable_modalities: Sequence[str],
    target_tissues: Sequence[str],
    head_mapping: Dict[str, Dict[str, Any]],
    gene_meta: Dict[str, Any],
    target_cell_types: Optional[Sequence[str]] = None,
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Scan candidate variants and compile AlphaGenome execution tasks for missing completable modalities.

    Returns:
        (triggers_list, warnings_list)
    """
    if target_cell_types is None:
        target_cell_types = ["K562"]
    triggers: List[Dict[str, Any]] = []
    warnings_out: List[str] = []
    counter = 0

    usable_contexts: Dict[str, List[str]] = {}
    for mod in completable_modalities:
        if mod == "eqtl":
            queried = _queried_contexts(evidence_by_key, "eqtl")
            usable: List[str] = []
            for t in target_tissues:
                if queried and t not in queried:
                    warnings_out.append(
                        f"CONFIG: eqtl target tissue '{t}' was not queried by Layer 1 "
                        f"(queried: {sorted(queried)}); no AlphaGenome jobs queued for it."
                    )
                else:
                    usable.append(t)
            usable_contexts[mod] = usable
        elif mod in ("mpra", "crispr"):
            usable_contexts[mod] = list(target_cell_types)
        else:
            usable_contexts[mod] = []

    fallback_reason = {
        "eqtl": "tissue_not_measured",
        "mpra": "cell_type_not_measured",
        "crispr": "cell_type_not_measured",
    }

    for v in variants:
        in_scope = filter_records_to_targets(
            evidence_by_key.get(v.key, []), target_tissues, target_cell_types
        )
        by_modality: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        for r in in_scope:
            by_modality[r["modality"]].append(r)

        for mod in completable_modalities:
            records = by_modality.get(mod, [])
            _, _, _, by_context = aggregate_modality_status(records)

            for ctx in usable_contexts.get(mod, []):
                # Context-specific record takes priority; otherwise fall back to context-less default
                rec = by_context.get(ctx, by_context.get("default"))
                if rec is not None and rec["flag"] != "absent":
                    continue  # Observed (positive or null) -> never complete
                if rec is not None and rec["category"] in NEVER_COMPLETE_CATEGORIES:
                    kind = "query failure" if rec["category"] == "ERROR_TRANSIENT" else "unrecognized reason"
                    warnings_out.append(
                        f"Variant {v.chrom}:{v.pos} {v.ref}>{v.alt} {mod} {ctx} had a "
                        f"{kind}: {rec.get('error') or rec.get('reason')}. "
                        f"Skipping AlphaGenome trigger."
                    )
                    continue

                counter += 1
                triggers.append(
                    _trigger_row(
                        counter,
                        v,
                        mod,
                        ctx,
                        gene_meta,
                        rec["reason"] if rec else fallback_reason.get(mod, "not_measured"),
                        rec["evidence_value"] if rec else {},
                        head_mapping=head_mapping,
                    )
                )

    return triggers, warnings_out


def compile_agreement_check_queue(
    variants: Sequence[Any],
    evidence_by_key: Dict[Tuple[str, int, str, str], List[Dict[str, Any]]],
    completable_modalities: Sequence[str],
    target_tissues: Optional[Sequence[str]] = None,
    target_cell_types: Optional[Sequence[str]] = None,
) -> List[Dict[str, Any]]:
    """Compile observed measurements on completable modalities into an agreement benchmark queue.

    Used strictly for model-vs-experiment benchmarking (never used in the fused score).
    """
    checks: List[Dict[str, Any]] = []
    counter = 0

    for v in variants:
        for r in filter_records_to_targets(
            evidence_by_key.get(v.key, []), target_tissues, target_cell_types
        ):
            mod, flag = r["modality"], r["missingness_flag"]
            if mod in completable_modalities and flag in ("observed_positive", "observed_null"):
                counter += 1
                checks.append(
                    {
                        "check_id": f"AG_AGREE_{counter:06d}",
                        "chrom": v.chrom,
                        "pos": v.pos,
                        "ref": v.ref,
                        "alt": v.alt,
                        "variant_key": f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}",
                        "modality": mod,
                        "context": _norm_ctx(r.get("context")),
                        "observed_flag": flag,
                        "observed_value": json.dumps(parse_evidence_value(r.get("evidence_value"))),
                        "validation_purpose": "internal_agreement_benchmark_only",
                    }
                )

    return checks
