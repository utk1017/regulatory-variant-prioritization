"""Data intake and synthetic fallback generation for Layer 2.

Discovers Layer 1 evidence outputs and Layer 0 candidate tables, loads ClinVar
labels, derives gene metadata, and provides a realistic synthetic fallback for
demonstration and pipeline regression tests when live inputs are not yet available.
"""
from __future__ import annotations

import json
import os
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd

from regvar.common.provenance import _file_hash
from regvar.layer1.variants import Variant, canon_key
from regvar.layer2.schema import validate_layer1_input


def to_variants(df: pd.DataFrame) -> List[Variant]:
    """Convert a DataFrame with chrom, pos, ref, alt into a deduplicated list of Variant namedtuples."""
    variants = []
    seen = set()
    for row in df[["chrom", "pos", "ref", "alt"]].itertuples(index=False):
        key = canon_key(row.chrom, row.pos, row.ref, row.alt)
        if key not in seen:
            seen.add(key)
            variants.append(Variant(key[0], key[1], key[2], key[3], key))
    return sorted(variants, key=lambda v: (v.chrom, v.pos, v.ref, v.alt))

SYNTH_ABSENT_REASON: Dict[str, str] = {
    "gwas": "not_in_catalog",
    "eqtl": "not_testable_rare_or_unknown_frequency",
    "gnomad_af": "not_observed_in_genome_callset",
    "mpra": "no_data_source_configured",
    "crispr": "no_data_source_configured",
}


def find_layer1_inputs(
    gene_symbol: str,
    explicit_csv: Optional[str] = None,
    explicit_manifest: Optional[str] = None,
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """Locate Layer 1 evidence CSV and manifest files across standard paths."""
    candidates = []
    if explicit_csv:
        manifest = explicit_manifest or os.path.join(os.path.dirname(explicit_csv), "layer1_manifest.json")
        candidates.append((explicit_csv, manifest, "explicit config"))
    candidates.append(
        (f"data/outputs/{gene_symbol}/layer1_evidence.csv", f"data/outputs/{gene_symbol}/layer1_manifest.json", "standard repo layout")
    )
    candidates.append(
        (f"layer1_outputs/{gene_symbol}_layer1_evidence.csv", None, "layer1 notebook layout")
    )

    for csv_path, manifest_path, label in candidates:
        if os.path.exists(csv_path):
            return csv_path, manifest_path, label
    return None, None, None


def find_layer0_inputs(
    gene_symbol: str, explicit_csv: Optional[str] = None
) -> Optional[str]:
    """Locate Layer 0 candidate truthset CSV across standard paths."""
    candidates = [
        f"data/outputs/{gene_symbol}/layer0_regulatory_truthset.csv",
        f"layer0_outputs/{gene_symbol}_layer0_regulatory_truthset.csv",
    ]
    if explicit_csv:
        candidates.insert(0, explicit_csv)
    for p in candidates:
        if os.path.exists(p):
            return p
    return None


def derive_fallback_gene_meta(variants: Sequence[Variant], gene_symbol: str) -> Dict[str, Any]:
    """Derive minimal gene_meta when no Layer 1 manifest is available."""
    if not variants:
        raise ValueError(
            f"Layer 1 evidence for {gene_symbol} loaded but produced zero candidate variants; "
            f"nothing to analyse. Check that Layer 1 ran successfully and its output isn't empty."
        )
    return {"gene_symbol": gene_symbol, "chrom": str(variants[0].chrom)}


def load_clinvar_labels(
    layer0_csv_path: Optional[str],
) -> Dict[Tuple[str, int, str, str], Dict[str, Any]]:
    """Load ClinVar clinical labels and review stars indexed by canonical variant key."""
    clinvar_labels: Dict[Tuple[str, int, str, str], Dict[str, Any]] = {}
    if layer0_csv_path and os.path.exists(layer0_csv_path):
        l0_df = pd.read_csv(layer0_csv_path, dtype={"chrom": str})
        for row in l0_df.to_dict("records"):
            k = canon_key(row["chrom"], row["pos"], row["ref"], row["alt"])
            clinvar_labels[k] = {
                "CLNSIG": row.get("CLNSIG"),
                "label": row.get("label"),
                "review_stars": row.get("review_stars"),
            }
    return clinvar_labels


def generate_synthetic_hbb_evidence() -> Tuple[
    pd.DataFrame,
    List[Variant],
    Dict[Tuple[str, int, str, str], Dict[str, Any]],
    Dict[str, Any],
    Dict[str, Any],
]:
    """Generate a realistic synthetic Layer 1 evidence table for HBB regulatory candidates."""
    synthetic_vars = [
        # 1. 5' UTR promoter variant: common benign (eQTL null, GWAS subthreshold, gnomAD AF 0.12)
        ("11", 5227002, "T", "C", "Benign", "benign", 2),
        # 2. Promoter pathogenic: rare beta-thalassemia variant (rare/unknown AF -> eQTL absent, GWAS absent)
        ("11", 5227172, "G", "A", "Pathogenic", "pathogenic", 2),
        # 3. LCR cCRE distal enhancer variant: tested in MPRA (observed_positive allelic skew)
        ("11", 5280000, "A", "G", "Pathogenic", "pathogenic", 1),
        # 4. Promoter variant with tested null MPRA and null eQTL
        ("11", 5227250, "C", "T", "Benign", "benign", 1),
        # 5. Deep non-coding candidate: zero coverage across all experimental modalities (Tier 3 sparse)
        ("11", 5229000, "G", "T", "Uncertain_significance", "vus", 0),
        # 6. Variant affected by transient query failure on GWAS
        ("11", 5226500, "A", "C", "Likely_benign", "benign", 1),
    ]

    variants = [Variant(c, p, r, a, canon_key(c, p, r, a)) for c, p, r, a, _, _, _ in synthetic_vars]
    clinvar_labels = {
        canon_key(c, p, r, a): {"CLNSIG": sig, "label": lbl, "review_stars": stars}
        for c, p, r, a, sig, lbl, stars in synthetic_vars
    }

    records: List[Dict[str, Any]] = []

    # Var 1: eQTL null in blood, GWAS null, gnomAD AF 0.12
    v1 = variants[0]
    records.append({
        "chrom": v1.chrom, "pos": v1.pos, "ref": v1.ref, "alt": v1.alt, "modality": "eqtl",
        "source_type": "experimental", "context": "Whole_Blood", "missingness_flag": "observed_null",
        "evidence_value": json.dumps({"reason": "common_variant_not_significant", "maf": 0.12}),
    })
    records.append({
        "chrom": v1.chrom, "pos": v1.pos, "ref": v1.ref, "alt": v1.alt, "modality": "eqtl",
        "source_type": "experimental", "context": "Liver", "missingness_flag": "observed_null",
        "evidence_value": json.dumps({"reason": "common_variant_not_significant", "maf": 0.12}),
    })
    records.append({
        "chrom": v1.chrom, "pos": v1.pos, "ref": v1.ref, "alt": v1.alt, "modality": "gwas",
        "source_type": "experimental", "context": "rs12345", "missingness_flag": "observed_null",
        "evidence_value": json.dumps({"reason": "subthreshold_association", "best_pvalue": 0.04}),
    })
    records.append({
        "chrom": v1.chrom, "pos": v1.pos, "ref": v1.ref, "alt": v1.alt, "modality": "gnomad_af",
        "source_type": "curated_clinical", "context": "genome", "missingness_flag": "observed_positive",
        "evidence_value": json.dumps({"af": 0.12}),
    })
    records.append({
        "chrom": v1.chrom, "pos": v1.pos, "ref": v1.ref, "alt": v1.alt, "modality": "mpra",
        "source_type": "experimental", "context": None, "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "not_tested_in_provided_table"}),
    })
    records.append({
        "chrom": v1.chrom, "pos": v1.pos, "ref": v1.ref, "alt": v1.alt, "modality": "crispr",
        "source_type": "experimental", "context": None, "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "no_data_source_configured"}),
    })

    # Var 2: Pathogenic promoter - rare, absent everywhere except ClinVar
    v2 = variants[1]
    records.append({
        "chrom": v2.chrom, "pos": v2.pos, "ref": v2.ref, "alt": v2.alt, "modality": "eqtl",
        "source_type": "experimental", "context": "Whole_Blood", "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "not_testable_rare_or_unknown_frequency", "maf": None}),
    })
    records.append({
        "chrom": v2.chrom, "pos": v2.pos, "ref": v2.ref, "alt": v2.alt, "modality": "eqtl",
        "source_type": "experimental", "context": "Liver", "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "not_testable_rare_or_unknown_frequency", "maf": None}),
    })
    records.append({
        "chrom": v2.chrom, "pos": v2.pos, "ref": v2.ref, "alt": v2.alt, "modality": "gwas",
        "source_type": "experimental", "context": None, "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "not_in_catalog"}),
    })
    records.append({
        "chrom": v2.chrom, "pos": v2.pos, "ref": v2.ref, "alt": v2.alt, "modality": "gnomad_af",
        "source_type": "curated_clinical", "context": "genome", "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "not_observed_in_genome_callset"}),
    })
    records.append({
        "chrom": v2.chrom, "pos": v2.pos, "ref": v2.ref, "alt": v2.alt, "modality": "mpra",
        "source_type": "experimental", "context": None, "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "not_tested_in_provided_table"}),
    })
    records.append({
        "chrom": v2.chrom, "pos": v2.pos, "ref": v2.ref, "alt": v2.alt, "modality": "crispr",
        "source_type": "experimental", "context": None, "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "no_data_source_configured"}),
    })

    # Var 3: Distal LCR enhancer variant - tested positive in MPRA
    v3 = variants[2]
    records.append({
        "chrom": v3.chrom, "pos": v3.pos, "ref": v3.ref, "alt": v3.alt, "modality": "mpra",
        "source_type": "experimental", "context": "K562", "missingness_flag": "observed_positive",
        "evidence_value": json.dumps({"log2_allelic_skew": 1.8, "pvalue": 1.2e-6, "dataset": "Tewhey2016"}),
    })
    records.append({
        "chrom": v3.chrom, "pos": v3.pos, "ref": v3.ref, "alt": v3.alt, "modality": "eqtl",
        "source_type": "experimental", "context": "Whole_Blood", "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "not_testable_rare_or_unknown_frequency", "maf": None}),
    })
    records.append({
        "chrom": v3.chrom, "pos": v3.pos, "ref": v3.ref, "alt": v3.alt, "modality": "gwas",
        "source_type": "experimental", "context": None, "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "not_in_catalog"}),
    })
    records.append({
        "chrom": v3.chrom, "pos": v3.pos, "ref": v3.ref, "alt": v3.alt, "modality": "gnomad_af",
        "source_type": "curated_clinical", "context": "genome", "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "not_observed_in_genome_callset"}),
    })
    records.append({
        "chrom": v3.chrom, "pos": v3.pos, "ref": v3.ref, "alt": v3.alt, "modality": "crispr",
        "source_type": "experimental", "context": None, "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "no_data_source_configured"}),
    })

    # Var 4: Promoter variant with tested null MPRA and null eQTL
    v4 = variants[3]
    records.append({
        "chrom": v4.chrom, "pos": v4.pos, "ref": v4.ref, "alt": v4.alt, "modality": "eqtl",
        "source_type": "experimental", "context": "Whole_Blood", "missingness_flag": "observed_null",
        "evidence_value": json.dumps({"reason": "common_variant_not_significant", "maf": 0.08}),
    })
    records.append({
        "chrom": v4.chrom, "pos": v4.pos, "ref": v4.ref, "alt": v4.alt, "modality": "gwas",
        "source_type": "experimental", "context": None, "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "not_in_catalog"}),
    })
    records.append({
        "chrom": v4.chrom, "pos": v4.pos, "ref": v4.ref, "alt": v4.alt, "modality": "gnomad_af",
        "source_type": "curated_clinical", "context": "genome", "missingness_flag": "observed_positive",
        "evidence_value": json.dumps({"af": 0.08}),
    })
    records.append({
        "chrom": v4.chrom, "pos": v4.pos, "ref": v4.ref, "alt": v4.alt, "modality": "mpra",
        "source_type": "experimental", "context": "K562", "missingness_flag": "observed_null",
        "evidence_value": json.dumps({"reason": "tested_not_significant", "pvalue": 0.42}),
    })
    records.append({
        "chrom": v4.chrom, "pos": v4.pos, "ref": v4.ref, "alt": v4.alt, "modality": "crispr",
        "source_type": "experimental", "context": None, "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "no_data_source_configured"}),
    })

    # Var 5: VUS - completely absent
    v5 = variants[4]
    for m in ("eqtl", "gwas", "gnomad_af", "mpra", "crispr"):
        records.append({
            "chrom": v5.chrom, "pos": v5.pos, "ref": v5.ref, "alt": v5.alt, "modality": m,
            "source_type": "experimental" if m != "gnomad_af" else "curated_clinical",
            "context": None, "missingness_flag": "absent",
            "evidence_value": json.dumps({"reason": SYNTH_ABSENT_REASON[m]}),
        })

    # Var 6: Transient query failure on GWAS
    v6 = variants[5]
    records.append({
        "chrom": v6.chrom, "pos": v6.pos, "ref": v6.ref, "alt": v6.alt, "modality": "gwas",
        "source_type": "experimental", "context": None, "missingness_flag": "absent",
        "evidence_value": json.dumps({"reason": "query_failed", "error": "HTTP 504 Gateway Timeout"}),
    })
    for m in ("eqtl", "gnomad_af", "mpra", "crispr"):
        records.append({
            "chrom": v6.chrom, "pos": v6.pos, "ref": v6.ref, "alt": v6.alt, "modality": m,
            "source_type": "experimental" if m != "gnomad_af" else "curated_clinical",
            "context": None, "missingness_flag": "absent",
            "evidence_value": json.dumps({"reason": SYNTH_ABSENT_REASON[m]}),
        })

    gene_meta = {
        "gene_symbol": "HBB", "gene_id": "ENSG00000244734", "chrom": "11",
        "start": 5225464, "end": 5229395, "strand": -1,
    }
    info = {
        "csv_path": "synthetic_in_memory", "csv_sha256": None,
        "layout": "synthetic fallback", "manifest_path": None,
    }

    return pd.DataFrame(records), variants, clinvar_labels, gene_meta, info
