"""Layer 2 Pipeline: Gap & Missingness Detection.

Orchestrates:
1. Intake of Layer 1 evidence matrix and Layer 0 candidate truth-set.
2. Disambiguation of absent evidence reasons (biological gaps vs. experimental nulls vs. transient errors).
3. Compilation of AlphaGenome computational trigger jobs for missing completable modalities (eQTL, MPRA, CRISPR).
4. Compilation of observed experimental hits into an agreement benchmark validation queue.
5. Assembly of per-variant Modality Availability vectors and calibration tiers (Tier 1, Tier 2, Tier 3).
6. Output serialization and write of provenance audit trail (``layer2_manifest.json``).

Outputs:
- ``data/outputs/<gene>/layer2_modality_availability.csv``
- ``data/outputs/<gene>/layer2_alphagenome_triggers.csv``
- ``data/outputs/<gene>/layer2_agreement_checks.csv``
- ``data/outputs/<gene>/layer2_manifest.json``
"""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import json
import os
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import pandas as pd

from regvar.common.config import load_config
from regvar.common.provenance import _file_hash, _git_commit_hash
from regvar.layer1.variants import canon_key
from regvar.layer2.intake import (
    derive_fallback_gene_meta,
    find_layer0_inputs,
    find_layer1_inputs,
    generate_synthetic_hbb_evidence,
    load_clinvar_labels,
    to_variants,
)
from regvar.layer2.reasons import (
    _norm_ctx,
    classify_evidence_reason,
    filter_records_to_targets,
    parse_evidence_value,
)
from regvar.layer2.schema import (
    LAYER2_AGREEMENT_COLUMNS,
    LAYER2_TRIGGER_COLUMNS,
    MODALITY_META,
    validate_layer1_input,
    validate_layer2_agreement_output,
    validate_layer2_availability_output,
    validate_layer2_triggers_output,
)
from regvar.layer2.tiers import (
    DEFAULT_TIER_THRESHOLDS,
    compute_evidence_tiers,
    label_coverage_warnings,
    tier_balance_warnings,
)
from regvar.layer2.triggers import (
    compile_agreement_check_queue,
    compile_alphagenome_triggers,
)

DEFAULT_L2: Dict[str, Any] = {
    "layer1_csv_path": None,
    "layer1_manifest_path": None,
    "layer0_csv_path": None,
    "target_tissues": ["Whole_Blood"],
    "target_cell_types": ["HEL92.1.7"],
    "completable_modalities": ["eqtl", "mpra", "crispr"],
    "regulatory_modalities": ["eqtl", "mpra", "crispr"],
    "tier_thresholds": DEFAULT_TIER_THRESHOLDS,
    "halt_on_query_failure": True,
    "strict": False,
}


def scan_query_errors(evidence_df: pd.DataFrame) -> List[str]:
    """Aggregate technical query failures across EVERY modality (including non-completable GWAS/gnomAD).
    A failed query is an unknown, never a biological absence."""
    from collections import defaultdict
    groups = defaultdict(lambda: defaultdict(int))
    for r in evidence_df.to_dict("records"):
        if r.get("missingness_flag") != "absent":
            continue
        ev = parse_evidence_value(r.get("evidence_value"))
        reason, err, errs = ev.get("reason"), ev.get("error"), ev.get("errors")
        if err or errs or classify_evidence_reason("absent", reason, None) == "ERROR_TRANSIENT":
            label = str(err or (errs[0] if errs else None) or reason)
            groups[r["modality"]][label] += 1
    out = []
    for mod in sorted(groups):
        detail = ", ".join(f"{k} x{n}" for k, n in sorted(groups[mod].items()))
        out.append(
            f"QUERY FAILURE: {mod}: {sum(groups[mod].values())} record(s) failed ({detail}). "
            f"These are unknowns, not biological absence: rerun Layer 1 with refresh_cache=True "
            f"and check the endpoint before trusting tiers or triggers."
        )
    return out


def enforce_data_quality_gate(blocking_warnings: List[str], halt: bool = True) -> str:
    """Return 'OK' or 'DEGRADED'; raise instead of silently producing degraded outputs when halt=True."""
    if not blocking_warnings:
        return "OK"
    if halt:
        raise RuntimeError(
            "Layer 2 halted: evidence queries failed, so tiers/triggers would be built on unknowns.\n  - "
            + "\n  - ".join(blocking_warnings)
            + "\nFix the failing source and rerun Layer 1, or set layer2.halt_on_query_failure: false in "
            "configs/hbb.yaml to continue with outputs explicitly marked DEGRADED."
        )
    return "DEGRADED"


def audit_context_scope(
    evidence_df: pd.DataFrame, target_tissues: Sequence[str], target_cell_types: Sequence[str]
) -> Tuple[List[str], pd.DataFrame]:
    """Report evidence that filter_records_to_targets() silently excludes (e.g. MPRA in K562 when the
    target is HEL92.1.7). Returns (warnings, per-modality context table)."""
    from collections import Counter
    warnings_out, rows = [], []
    for mod, targets in (("eqtl", target_tissues), ("mpra", target_cell_types), ("crispr", target_cell_types)):
        sub = evidence_df[evidence_df["modality"] == mod]
        ctxs = Counter(_norm_ctx(c) for c in sub["context"])
        named = {c: n for c, n in ctxs.items() if c != "default"}
        inside = {c: n for c, n in named.items() if c in targets}
        outside = {c: n for c, n in named.items() if c not in targets}
        rows.append({
            "modality": mod,
            "targets": ";".join(targets),
            "contexts_in_evidence": ";".join(f"{c}:{n}" for c, n in sorted(named.items())) or "(none)",
            "records_in_scope": int(sum(inside.values())),
            "records_excluded": int(sum(outside.values())),
        })
        if outside:
            warnings_out.append(
                f"SCOPE: {mod} evidence in context(s) {sorted(outside)} is outside the configured target(s) "
                f"{list(targets)} and is EXCLUDED from tiers and triggers. If that is the wrong target, "
                f"fix layer2.target_{'tissues' if mod == 'eqtl' else 'cell_types'} in configs/hbb.yaml."
            )
        if mod in ("mpra", "crispr") and named and not inside:
            warnings_out.append(
                f"SCOPE: NONE of the {mod} table's cell types match the target "
                f"{list(targets)} -- all {mod} evidence is being ignored."
            )
    return warnings_out, pd.DataFrame(rows)


def tier1_feasibility(availability_df: pd.DataFrame, min_modalities: int = 2) -> pd.DataFrame:
    """Explain Tier 1 sparsity and list the variants a second allele-resolved source would promote."""
    allele_cols = [c for c in ("eqtl_status", "mpra_status", "crispr_status") if c in availability_df.columns]
    observed = availability_df[allele_cols].isin(["observed_positive", "observed_null"])
    n_obs = observed.sum(axis=1)
    cand = availability_df[(n_obs == min_modalities - 1) & (n_obs >= 1)].copy()
    if not cand.empty:
        cand["has_modalities"] = observed.loc[cand.index].apply(
            lambda r: ";".join(c.replace("_status", "") for c, v in r.items() if v), axis=1
        )
        cand["needs_additional_modalities"] = 1
    else:
        cand["has_modalities"] = pd.Series(dtype=str)
        cand["needs_additional_modalities"] = pd.Series(dtype=int)
    return cand


def resolve_out_dir(gene_symbol: str, is_synthetic: bool) -> str:
    """Real runs write to data/outputs/<gene>/; synthetic demos are segregated to synthetic_demo/."""
    base = os.path.join("data", "outputs", gene_symbol)
    return os.path.join(base, "synthetic_demo") if is_synthetic else base


def write_layer2_manifest(out_dir: str, payload: Dict[str, Any]) -> str:
    """Serialize Layer 2 run provenance and audit trail to JSON."""
    os.makedirs(out_dir, exist_ok=True)
    manifest_path = os.path.join(out_dir, "layer2_manifest.json")
    payload["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    payload["git_commit"] = _git_commit_hash()
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, default=str)
    return manifest_path


def run_layer2(
    config_input: Union[str, Dict[str, Any]]
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Execute Layer 2 Gap Detection and AlphaGenome Trigger Compilation.

    Args:
        config_input: Path to config YAML file or pre-loaded config dict.

    Returns:
        (availability_df, triggers_df, agreement_df)
    """
    if isinstance(config_input, str):
        config = load_config(config_input) if os.path.exists(config_input) else {}
        config_path = config_input
    else:
        config = config_input or {}
        config_path = None

    gene_symbol = config.get("gene_symbol", "HBB")
    l1_cfg = config.get("layer1") or {}
    l2_cfg_user = config.get("layer2") or {}

    # Merge effective Layer 2 configuration
    l2_cfg = dict(DEFAULT_L2)
    l2_cfg.update(l2_cfg_user)

    # If target_tissues not explicitly set in layer2, default to layer1.gtex_tissues
    if "target_tissues" not in l2_cfg_user and "gtex_tissues" in l1_cfg:
        l2_cfg["target_tissues"] = list(l1_cfg["gtex_tissues"])

    # If paths provided at top level, preserve them
    if config.get("layer1_csv_path"):
        l2_cfg["layer1_csv_path"] = config["layer1_csv_path"]
    if config.get("layer1_manifest_path"):
        l2_cfg["layer1_manifest_path"] = config["layer1_manifest_path"]
    if config.get("layer0_csv_path"):
        l2_cfg["layer0_csv_path"] = config["layer0_csv_path"]
    if "strict" in config:
        l2_cfg["strict"] = config["strict"]

    print(f"Running Layer 2 Gap Detection & Trigger Compilation for {gene_symbol}...")

    # Discover inputs
    l1_csv, l1_manifest, layout = find_layer1_inputs(
        gene_symbol, l2_cfg.get("layer1_csv_path"), l2_cfg.get("layer1_manifest_path")
    )
    l0_csv = find_layer0_inputs(gene_symbol, l2_cfg.get("layer0_csv_path"))

    clinvar_labels = load_clinvar_labels(l0_csv)
    if clinvar_labels:
        print(f"Loaded {len(clinvar_labels)} variant annotations from Layer 0 ({l0_csv})")

    if l1_csv:
        print(f"Loading live Layer 1 evidence from {l1_csv} ({layout})...")
        evidence_df = pd.read_csv(l1_csv, dtype={"chrom": str})
        validate_layer1_input(evidence_df)
        is_synthetic = False
        variants = to_variants(evidence_df)

        gene_meta = derive_fallback_gene_meta(variants, gene_symbol)
        if l1_manifest and os.path.exists(l1_manifest):
            with open(l1_manifest, "r", encoding="utf-8") as f:
                m = json.load(f)
            gene_meta = m.get("gene_metadata", gene_meta)
        layer1_info = {
            "csv_path": l1_csv,
            "csv_sha256": _file_hash(l1_csv),
            "layout": layout,
            "manifest_path": l1_manifest,
        }
        print(f"Loaded {len(evidence_df)} evidence records across {len(variants)} candidates.")
    else:
        if l2_cfg.get("strict"):
            raise FileNotFoundError(f"Strict mode enabled: Layer 1 evidence CSV not found for {gene_symbol}.")
        print("INFO: Live Layer 1 outputs not detected on disk yet.")
        print("!" * 70)
        print("!! SYNTHETIC DEMO MODE: live Layer 1 inputs missing. Generating synthetic demonstration.")
        print("!! Outputs will be segregated to synthetic_demo/ and flagged is_synthetic=true.")
        print("!" * 70)
        is_synthetic = True
        evidence_df, variants, synth_clinvar, gene_meta, layer1_info = generate_synthetic_hbb_evidence()
        if not clinvar_labels:
            clinvar_labels = synth_clinvar
        print(f"Loaded {len(evidence_df)} records for {len(variants)} synthetic {gene_symbol} candidate variants.")

    # Index evidence records by canonical key
    evidence_by_key: Dict[Tuple[str, int, str, str], List[Dict[str, Any]]] = defaultdict(list)
    for rec in evidence_df.to_dict("records"):
        k = canon_key(rec["chrom"], rec["pos"], rec["ref"], rec["alt"])
        evidence_by_key[k].append(rec)

    # 1. Modality Availability & Stratification
    availability_df = compute_evidence_tiers(
        variants=variants,
        evidence_by_key=evidence_by_key,
        clinvar_labels=clinvar_labels,
        regulatory_modalities=l2_cfg.get("regulatory_modalities", ["eqtl", "mpra", "crispr"]),
        target_tissues=l2_cfg["target_tissues"],
        target_cell_types=l2_cfg["target_cell_types"],
        tier_thresholds=l2_cfg.get("tier_thresholds"),
    )
    validate_layer2_availability_output(availability_df)

    # 2. AlphaGenome Trigger Queue
    completable = l2_cfg.get("completable_modalities") or [
        m for m, meta in MODALITY_META.items() if meta.get("completable")
    ]
    triggers, transient_warnings = compile_alphagenome_triggers(
        variants=variants,
        evidence_by_key=evidence_by_key,
        completable_modalities=completable,
        target_tissues=l2_cfg["target_tissues"],
        head_mapping=MODALITY_META,
        gene_meta=gene_meta,
        target_cell_types=l2_cfg["target_cell_types"],
    )
    triggers_df = pd.DataFrame(triggers, columns=LAYER2_TRIGGER_COLUMNS)
    validate_layer2_triggers_output(triggers_df, completable_modalities=completable)

    # 3. Agreement Check Queue
    agreements = compile_agreement_check_queue(
        variants=variants,
        evidence_by_key=evidence_by_key,
        completable_modalities=completable,
        target_tissues=l2_cfg["target_tissues"],
        target_cell_types=l2_cfg["target_cell_types"],
    )
    agreement_df = pd.DataFrame(agreements, columns=LAYER2_AGREEMENT_COLUMNS)
    validate_layer2_agreement_output(agreement_df)

    tier_warnings = tier_balance_warnings(availability_df)
    label_warnings = label_coverage_warnings(availability_df)

    # 4. Diagnostics & Data Quality Gate
    blocking_warnings = scan_query_errors(evidence_df)
    scope_warnings, scope_table = audit_context_scope(evidence_df, l2_cfg["target_tissues"], l2_cfg["target_cell_types"])
    data_quality = enforce_data_quality_gate(
        blocking_warnings, halt=l2_cfg.get("halt_on_query_failure", True) and not is_synthetic
    )

    tier_1_min = l2_cfg.get("tier_thresholds", {}).get("tier_1_min_modalities", 2)
    tier1_candidates_df = tier1_feasibility(availability_df, min_modalities=tier_1_min)

    print(f"\nLayer 2 Assembly Complete:")
    print(f" - Modality Availability Matrix: {len(availability_df)} variants")
    print(f" - AlphaGenome Trigger Queue:    {len(triggers_df)} inference jobs")
    print(f" - Agreement Check Queue:        {len(agreement_df)} validation benchmarks")
    print(f" - Pipeline Warnings (query/config): {len(transient_warnings)}")
    print(f" - Tier Balance Warnings:        {len(tier_warnings)}")
    print(f" - Label Coverage Warnings:      {len(label_warnings)}")
    print(f" - Scope Warnings:               {len(scope_warnings)}")
    print(f" - Data Quality Gate:            {data_quality}")

    # Save outputs
    out_dir = resolve_out_dir(gene_symbol, is_synthetic)
    os.makedirs(out_dir, exist_ok=True)

    avail_path = os.path.join(out_dir, "layer2_modality_availability.csv")
    triggers_path = os.path.join(out_dir, "layer2_alphagenome_triggers.csv")
    agreement_path = os.path.join(out_dir, "layer2_agreement_checks.csv")
    tier1_path = os.path.join(out_dir, "layer2_tier1_candidates_needing_second_source.csv")

    availability_df.to_csv(avail_path, index=False)
    triggers_df.to_csv(triggers_path, index=False)
    agreement_df.to_csv(agreement_path, index=False)
    tier1_candidates_df.to_csv(tier1_path, index=False)

    tier_counts = {k: int(v) for k, v in availability_df["evidence_tier"].value_counts().items()}
    trigger_by_modality = (
        {k: int(v) for k, v in triggers_df["modality"].value_counts().items()}
        if not triggers_df.empty
        else {}
    )
    trigger_by_context = (
        {f"{m}|{c}": int(n) for (m, c), n in triggers_df.groupby(["modality", "context"]).size().items()}
        if not triggers_df.empty
        else {}
    )
    trigger_by_output_type: Dict[str, int] = {}
    if not triggers_df.empty:
        for s in triggers_df["alphagenome_output_types"]:
            for t in str(s).split(";"):
                trigger_by_output_type[t] = trigger_by_output_type.get(t, 0) + 1

    manifest_payload: Dict[str, Any] = {
        "layer": 2,
        "is_synthetic": bool(is_synthetic),
        "data_quality": data_quality,
        "gene_symbol": gene_symbol,
        "gene_metadata": gene_meta,
        "layer1_input": layer1_info,
        "config_path": config_path,
        "config_sha256": _file_hash(config_path) if config_path else None,
        "layer2_effective_config": l2_cfg,
        "modality_meta": MODALITY_META,
        "n_candidate_variants": len(variants),
        "tier_distribution": tier_counts,
        "tier1_feasibility": {
            "tier1_count": int((availability_df["evidence_tier"] == "Tier 1 (Evidence-Rich)").sum()),
            "promotable_candidates_needing_second_source": len(tier1_candidates_df),
        },
        "alphagenome_triggers": {
            "total_jobs": len(triggers_df),
            "by_modality": trigger_by_modality,
            "by_context": trigger_by_context,
            "by_output_type": trigger_by_output_type,
            "binding_status": "unverified: Layer 3 must map output types to real AlphaGenome API enums",
        },
        "agreement_checks": {
            "total_checks": len(agreement_df),
        },
        "pipeline_warnings_count": len(transient_warnings),
        "pipeline_warnings": transient_warnings,
        "scope_warnings": scope_warnings,
        "tier_warnings": tier_warnings,
        "label_warnings": label_warnings,
    }

    manifest_path = write_layer2_manifest(out_dir, manifest_payload)

    print(f"Saved Modality Availability to: {avail_path}")
    print(f"Saved AlphaGenome Triggers to:  {triggers_path}")
    print(f"Saved Agreement Checks to:      {agreement_path}")
    print(f"Saved Tier 1 Candidates to:     {tier1_path}")
    print(f"Saved Manifest to:              {manifest_path}")

    return availability_df, triggers_df, agreement_df
