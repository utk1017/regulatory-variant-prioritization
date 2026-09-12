import os
import pandas as pd
from regvar.common.config import load_config
from regvar.layer0.gene_meta import fetch_gene_metadata
from regvar.layer0.reference import download_reference
from regvar.layer0.clinvar import (
    discover_latest_clinvar_date, download_clinvar_snapshot, 
    extract_region, normalize_vcf, extract_clinvar_labels
)
from regvar.layer0.intervals import fetch_ccres, flag_interval_overlap
from regvar.layer0.vep import vep_annotate_batch, parse_vep_results

REGULATORY_UTR_TERMS = [
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "regulatory_region_variant",
    "TF_binding_site_variant",
]

def build_evidence(row):
    evidence = [f"interval:{label}" for label in row.get("interval_evidence", [])]
    
    terms = row.get("vep_all_transcript_terms", [])
    if not isinstance(terms, list): terms = []
    reg_terms = row.get("vep_regulatory_terms", [])
    if not isinstance(reg_terms, list): reg_terms = []
    
    utr_terms = [t for t in (terms + reg_terms) if any(k in t for k in REGULATORY_UTR_TERMS)]
    evidence += [f"vep:{term}" for term in utr_terms]
    return evidence

def run_layer0(config_path):
    config = load_config(config_path)
    
    gene_symbol = config["gene_symbol"]
    canonical_transcript = config["canonical_transcript"]
    flank_bp = config.get("flank_bp", 100000)
    promoter_window_bp = config.get("promoter_window_bp", 2500)
    clinvar_date = config.get("clinvar_snapshot_date")
    
    print(f"Running Layer 0 for {gene_symbol}...")
    
    gene_meta = fetch_gene_metadata(gene_symbol)
    
    if gene_meta["strand"] == 1:
        promoter_start = gene_meta["start"] - promoter_window_bp
        promoter_end = gene_meta["start"]
    else:
        promoter_start = gene_meta["end"]
        promoter_end = gene_meta["end"] + promoter_window_bp
        
    flank_start = max(1, gene_meta["start"] - flank_bp)
    flank_end = gene_meta["end"] + flank_bp
    region_str = f"{gene_meta['chrom']}:{flank_start}-{flank_end}"
    
    ref_fasta = download_reference(gene_meta["chrom"])
    
    if not clinvar_date:
        clinvar_date = discover_latest_clinvar_date()
        print(f"Auto-discovered ClinVar date: {clinvar_date}")
        
    clinvar_full = download_clinvar_snapshot(clinvar_date)
    clinvar_region = extract_region(clinvar_full, region_str, "data/cache/clinvar_region.vcf.gz")
    clinvar_norm = normalize_vcf(clinvar_region, ref_fasta, "data/cache/clinvar_normalized.vcf.gz")
    
    variants_df = extract_clinvar_labels(clinvar_norm)
    if variants_df.empty:
        print("No variants found passing filters in this region.")
        return variants_df
        
    ccre_df = fetch_ccres(gene_meta["chrom"], flank_start, flank_end)
    promoter_interval_df = pd.DataFrame({
        "Chromosome": [gene_meta["chrom"]],
        "Start": [promoter_start],
        "End": [promoter_end],
        "ccre_label": ["computed_promoter_window"],
    })
    reg_intervals_df = pd.concat([ccre_df, promoter_interval_df], ignore_index=True)
    
    variants_df = flag_interval_overlap(variants_df, reg_intervals_df)
    
    vep_results = vep_annotate_batch(variants_df, 
                                     batch_size=config.get("vep_batch_size", 200), 
                                     timeout=config.get("vep_timeout_s", 30), 
                                     max_retries=config.get("vep_max_retries", 3))
    vep_df = parse_vep_results(vep_results, canonical_transcript, gene_meta["gene_id"])
    
    df_merged = variants_df.merge(vep_df, on=["chrom", "pos", "ref", "alt"], how="left")
    df_merged["regulatory_evidence"] = df_merged.apply(build_evidence, axis=1)
    df_merged["is_regulatory"] = df_merged["regulatory_evidence"].str.len() > 0
    
    final_df = df_merged[df_merged["is_regulatory"]].drop_duplicates(subset=["chrom", "pos", "ref", "alt"]).reset_index(drop=True)
    
    out_dir = "data/outputs"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{gene_symbol}_layer0_regulatory_truthset.csv")
    final_df.to_csv(out_path, index=False)
    print(f"Layer 0 complete. Saved {len(final_df)} regulatory variants to {out_path}")
    
    return final_df
