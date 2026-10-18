"""Layer 1 Pipeline: Multi-source evidence collection and harmonization.

Orchestrates evidence extraction across:
1. gnomAD (population AF & MAF proxy)
2. GWAS Catalog (windowed lead SNP associations)
3. GTEx Portal (eQTL expression effects by tissue)
4. MPRA (reporter assay screens from local table)
5. CRISPR (perturbation screens from local table)

Outputs:
- ``data/outputs/<gene>/layer1_evidence.csv``: complete matrix of EvidenceRecords.
- ``data/outputs/<gene>/layer1_manifest.json``: run provenance and audit trail.
"""
from __future__ import annotations

import json
import os
from typing import Any, Dict, List, Optional

import pandas as pd

from regvar.common.config import load_config
from regvar.common.provenance import _file_hash, _git_commit_hash
from regvar.layer1.cache import HTTPCache, make_fetch
from regvar.layer1.gnomad import (
    build_gnomad_records,
    collect_gnomad_evidence,
    maf_from_gnomad,
)
from regvar.layer1.gtex import build_eqtl_records, collect_eqtl_evidence
from regvar.layer1.gwas import build_gwas_records, collect_gwas_evidence
from regvar.layer1.intake import load_layer0_output
from regvar.layer1.local_tables import (
    CRISPR_COLUMNS,
    MPRA_COLUMNS,
    build_local_table_records,
    crispr_value,
    load_local_evidence_table,
    mpra_value,
)
from regvar.layer1.schema import (
    MODALITIES,
    MODALITY_META,
    records_to_dataframe,
    validate_layer1_output,
    validate_matrix_complete,
)
from regvar.layer1.variants import to_variants


DEFAULT_L1 = {
    "layer0_csv_path": None,
    "strict": False,
    "use_http_cache": True,
    "http_cache_dir": "data/cache/http",
    "cache_ttl_days": 7,
    "refresh_cache": False,
    "gtex_dataset": "gtex_v8",
    "gtex_gencode_version": "v26",
    "gtex_tissues": ["Whole_Blood"],
    "gtex_items_per_page": 250,
    "gtex_max_pages": 40,
    "eqtl_testable_min_maf": 0.01,
    "gwas_p_threshold": 5.0e-8,
    "gwas_window_bp": 25000,
    "gnomad_dataset": "gnomad_r4",
    "gnomad_window_bp": 25000,
    "gnomad_request_delay_s": 6.0,
    "mpra_evidence_path": None,
    "crispr_evidence_path": None,
    "sig_pvalue_threshold": 0.05,
}


def write_layer1_manifest(out_dir: str, payload: dict) -> str:
    from datetime import datetime, timezone

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit_hash(),
    }
    manifest.update(payload)
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "layer1_manifest.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, default=str)
    return path


def run_layer1(config_path: str, fetch: Optional[Any] = None) -> Optional[pd.DataFrame]:
    """Execute Layer 1 multi-source evidence collection for a gene configuration."""
    config = load_config(config_path) if os.path.exists(config_path) else {}
    gene_symbol = config.get("gene_symbol", "HBB")
    user_l1 = config.get("layer1") or {}
    l1_cfg = {**DEFAULT_L1, **user_l1}

    print(f"Running Layer 1 evidence collection for {gene_symbol}...")

    # Initialize cache and fetch closure if not provided
    cache = HTTPCache(
        cache_dir=l1_cfg["http_cache_dir"],
        enabled=l1_cfg["use_http_cache"],
        ttl_days=l1_cfg.get("cache_ttl_days", 7),
        refresh=l1_cfg.get("refresh_cache", False),
    )
    if fetch is None:
        fetch = make_fetch(cache)

    try:
        variants_df, gene_meta, layer0_info = load_layer0_output(
            gene_symbol, l1_cfg["layer0_csv_path"], fetch=fetch
        )
    except FileNotFoundError as err:
        if l1_cfg.get("strict"):
            raise
        print(f"INFO: Layer 0 outputs not found yet ({err}). Skipping live collection.")
        return None

    variants = to_variants(variants_df)
    if not variants:
        print("Layer 0 produced 0 candidate variants. Nothing to collect.")
        return pd.DataFrame()

    chrom = str(gene_meta["chrom"]).replace("chr", "")
    if {v.chrom for v in variants} != {chrom}:
        raise ValueError(
            "Candidate variants span more than the gene's chromosome; "
            "multi-chromosome candidate sets are not supported."
        )

    run_warnings: List[str] = []

    # 1. gnomAD
    print("Collecting gnomAD allele frequencies...")
    gnomad = collect_gnomad_evidence(
        variants,
        l1_cfg["gnomad_dataset"],
        l1_cfg["gnomad_window_bp"],
        l1_cfg["gnomad_request_delay_s"],
        fetch=fetch,
    )
    maf_by_key = maf_from_gnomad(gnomad)
    for w in gnomad["windows"]:
        if w["error"]:
            run_warnings.append(
                f"gnomAD window {w['chrom']}:{w['start']}-{w['stop']} failed: {w['error']}"
            )
    print(
        f"gnomAD: {len(gnomad['windows'])} window(s), "
        f"{len(gnomad['by_key'])} variant(s) with genome AF"
    )

    # 2. GWAS Catalog (Windowed)
    print("Collecting GWAS Catalog associations...")
    positions = sorted({v.pos for v in variants})
    gwas = collect_gwas_evidence(
        chrom,
        positions,
        fetch=fetch,
        window_bp=l1_cfg.get("gwas_window_bp", 25000),
    )
    failed_gwas_wins = [w for w in gwas["windows"] if w["error"]]
    if failed_gwas_wins:
        run_warnings.append(
            f"GWAS: {len(failed_gwas_wins)} of {len(gwas['windows'])} window(s) failed query"
        )
    print(
        f"GWAS: {gwas['n_snps_in_region']} catalog SNP(s) across {len(gwas['windows'])} window(s), "
        f"{len(gwas['by_pos'])} at candidate positions"
    )

    # 3. GTEx eQTL
    print("Collecting GTEx eQTLs...")
    tissues = list(l1_cfg["gtex_tissues"])
    eqtl = collect_eqtl_evidence(
        gene_meta,
        tissues,
        l1_cfg["gtex_dataset"],
        l1_cfg["gtex_gencode_version"],
        l1_cfg["gtex_items_per_page"],
        l1_cfg["gtex_max_pages"],
        fetch=fetch,
    )
    if eqtl["resolve_error"]:
        run_warnings.append(f"GTEx gene resolution failed: {eqtl['resolve_error']}")
    if eqtl["resolve_note"]:
        run_warnings.append(f"GTEx gene resolution note: {eqtl['resolve_note']}")
    for t, info in eqtl["tissues"].items():
        if info["error"]:
            run_warnings.append(f"GTEx tissue {t}: {info['error']} after {info['n_rows']} row(s)")
        elif info["n_rows"] == 0:
            run_warnings.append(
                f"GTEx tissue {t}: 0 significant eQTL rows returned for {gene_meta['gene_symbol']}"
            )
    print(
        f"GTEx: gencode_id={eqtl['gencode_id']}, "
        + ", ".join(f"{t}={i['n_rows']} rows" for t, i in eqtl["tissues"].items())
    )

    # 4. Local tables (MPRA & CRISPR)
    mpra_index = load_local_evidence_table(
        l1_cfg["mpra_evidence_path"], MPRA_COLUMNS, "MPRA", ("log2_allelic_skew",)
    )
    crispr_index = load_local_evidence_table(
        l1_cfg["crispr_evidence_path"], CRISPR_COLUMNS, "CRISPR", ("effect_size",)
    )

    # 5. Assemble all records
    sig_p = l1_cfg["sig_pvalue_threshold"]
    all_records = (
        build_gwas_records(gwas, variants, l1_cfg["gwas_p_threshold"])
        + build_eqtl_records(eqtl, variants, tissues, maf_by_key, l1_cfg["eqtl_testable_min_maf"])
        + build_gnomad_records(gnomad, variants, l1_cfg["gnomad_dataset"])
        + build_local_table_records("mpra", variants, mpra_index, sig_p, mpra_value)
        + build_local_table_records("crispr", variants, crispr_index, sig_p, crispr_value)
    )

    evidence_df = records_to_dataframe(all_records)
    validate_layer1_output(evidence_df)
    validate_matrix_complete(evidence_df, variants)

    # 6. Save outputs & manifest
    out_dir = os.path.join("data/outputs", gene_symbol)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "layer1_evidence.csv")
    evidence_df.to_csv(out_path, index=False)

    flag_counts = {
        m: {f: int(n) for f, n in g["missingness_flag"].value_counts().items()}
        for m, g in evidence_df.groupby("modality")
    }

    manifest_path = write_layer1_manifest(
        out_dir,
        {
            "layer": 1,
            "gene_symbol": gene_symbol,
            "gene_metadata": gene_meta,
            "layer0_input": layer0_info,
            "config_path": config_path,
            "config_sha256": _file_hash(config_path),
            "layer1_effective_config": l1_cfg,
            "modality_meta": MODALITY_META,
            "n_candidate_variants": len(variants),
            "n_evidence_records": len(evidence_df),
            "flag_counts": flag_counts,
            "cache": {
                "enabled": l1_cfg["use_http_cache"],
                "cache_dir": l1_cfg["http_cache_dir"],
                "hits": cache.hits,
                "misses": cache.misses,
                "oldest_entry_utc": cache.oldest_entry_utc,
            },
            "gtex": {
                "gencode_id": eqtl["gencode_id"],
                "dataset": l1_cfg["gtex_dataset"],
                "tissues": {
                    t: {"n_rows": i["n_rows"], "error": i["error"]}
                    for t, i in eqtl["tissues"].items()
                },
            },
            "gnomad": {"dataset": l1_cfg["gnomad_dataset"], "windows": gnomad["windows"]},
            "gwas": {"windows": gwas["windows"], "n_snps_in_span": gwas["n_snps_in_region"]},
            "mpra_evidence_path": l1_cfg["mpra_evidence_path"],
            "mpra_sha256": _file_hash(l1_cfg["mpra_evidence_path"]),
            "crispr_evidence_path": l1_cfg["crispr_evidence_path"],
            "crispr_sha256": _file_hash(l1_cfg["crispr_evidence_path"]),
            "any_query_failed": bool(run_warnings),
            "warnings": run_warnings,
        },
    )

    print(f"Layer 1 complete. Saved {len(evidence_df)} evidence records to {out_path}")
    print(f"Provenance manifest: {manifest_path}")
    if run_warnings:
        print(f"{len(run_warnings)} warning(s) occurred:")
        for w in run_warnings:
            print("  -", w)

    return evidence_df
