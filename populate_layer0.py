import os

base = r"d:\projects\regulatory-variant-prioritization"

files = {
    "src/regvar/common/config.py": """import yaml

def load_config(config_path):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)
""",

    "src/regvar/common/utils.py": """import subprocess

def _run(cmd, check=True):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        msg = f"Command failed (exit {r.returncode}): {cmd}\\nSTDERR:\\n{r.stderr}"
        if check:
            raise RuntimeError(msg)
        else:
            print(msg)
    return r

def to_ucsc_chrom(chrom):
    chrom = str(chrom)
    return chrom if chrom.startswith("chr") else f"chr{chrom}"

def to_plain_chrom(chrom):
    chrom = str(chrom)
    return chrom[3:] if chrom.startswith("chr") else chrom
""",

    "src/regvar/layer0/gene_meta.py": """import requests

def fetch_gene_metadata(gene_symbol, species="homo_sapiens", timeout=15):
    url = f"https://rest.ensembl.org/lookup/symbol/{species}/{gene_symbol}"
    r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=timeout)
    r.raise_for_status()
    data = r.json()
    return {
        "gene_symbol": gene_symbol,
        "gene_id": data["id"],
        "chrom": str(data["seq_region_name"]),
        "start": int(data["start"]),
        "end": int(data["end"]),
        "strand": int(data["strand"]),
    }
""",

    "src/regvar/layer0/reference.py": """import os
from regvar.common.utils import _run

def download_reference(chrom, out_dir="data/cache", out_prefix="ref"):
    os.makedirs(out_dir, exist_ok=True)
    fa_gz = os.path.join(out_dir, f"{out_prefix}_chr{chrom}.fa.gz")
    fa = os.path.join(out_dir, f"{out_prefix}_chr{chrom}.fa")
    url = f"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr{chrom}.fa.gz"

    _run(f"wget -nc {url} -O {fa_gz}")
    integrity = _run(f"gzip -t {fa_gz}", check=False)
    if integrity.returncode != 0:
        print(f"Corrupt download detected, re-fetching {fa_gz}")
        _run(f"rm -f {fa_gz}")
        _run(f"wget {url} -O {fa_gz}")
        _run(f"gzip -t {fa_gz}")

    _run(f"gunzip -kf {fa_gz}")
    _run(f"sed -i 's/^>chr{chrom}/>{chrom}/' {fa}")
    _run(f"samtools faidx {fa}")
    if not os.path.exists(f"{fa}.fai"):
        raise RuntimeError(f"samtools faidx did not produce an index for {fa}")
    return fa
""",

    "src/regvar/layer0/clinvar.py": """import os
import re
import requests
import pysam
import pandas as pd
from regvar.common.utils import _run

REVIEW_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "no_assertion_criteria_provided": 0,
    "no_assertion_provided": 0,
}

def discover_latest_clinvar_date(timeout=15):
    r = requests.get("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/", timeout=timeout)
    dates = sorted(set(re.findall(r'clinvar_(\d{8})\.vcf\.gz"', r.text)), reverse=True)
    if not dates:
        raise RuntimeError("Could not find any dated ClinVar snapshot")
    return dates[0]

def clinvar_archive_url(date_str):
    return f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_{date_str}.vcf.gz"

def download_clinvar_snapshot(date_str, out_name="data/cache/clinvar_snapshot.vcf.gz"):
    os.makedirs(os.path.dirname(out_name), exist_ok=True)
    url = clinvar_archive_url(date_str)
    _run(f"rm -f {out_name}")
    _run(f"wget {url} -O {out_name}")
    _run(f"gzip -t {out_name}")
    _run(f"tabix -f -p vcf {out_name}")
    return out_name

def extract_region(vcf_gz, region_str, out_name):
    _run(f"bcftools view {vcf_gz} -r {region_str} -Oz -o {out_name}")
    _run(f"tabix -f -p vcf {out_name}")
    return out_name

def normalize_vcf(input_vcf_gz, ref_fasta, out_name):
    split_name = "data/cache/split_tmp.vcf.gz"
    _run(f"bcftools norm -f {ref_fasta} -m -both {input_vcf_gz} -Oz -o {split_name}")
    _run(f"tabix -f -p vcf {split_name}")
    _run(f"bcftools norm -f {ref_fasta} -d exact {split_name} -Oz -o {out_name}")
    _run(f"tabix -f -p vcf {out_name}")
    return out_name

def extract_clinvar_labels(vcf_path):
    vf = pysam.VariantFile(vcf_path)
    rows = []
    for rec in vf:
        clnsig = rec.info.get("CLNSIG", ("",))
        clnsig = ",".join(clnsig) if isinstance(clnsig, tuple) else str(clnsig)
        if "Conflicting" in clnsig:
            continue
        if "Pathogenic" in clnsig:
            label = "pathogenic"
        elif "Benign" in clnsig:
            label = "benign"
        else:
            continue

        revstat = rec.info.get("CLNREVSTAT", ("",))
        revstat = ",".join(revstat) if isinstance(revstat, tuple) else str(revstat)
        stars = max([v for k, v in REVIEW_STARS.items() if k in revstat], default=0)

        for alt in (rec.alts or []):
            rows.append({
                "chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref, "alt": alt,
                "CLNSIG": clnsig, "CLNREVSTAT": revstat, "review_stars": stars,
                "label": label,
            })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt"]).reset_index(drop=True)
        df = df[df["review_stars"] >= 1].reset_index(drop=True)
    return df
""",

    "src/regvar/layer0/intervals.py": """import os
import pandas as pd
import pyranges as pr
from regvar.common.utils import _run, to_ucsc_chrom, to_plain_chrom

def fetch_bigbedtobed_tool():
    tool_path = "scripts/bigBedToBed"
    if not os.path.exists(tool_path):
        _run(f"wget -q http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed -O {tool_path}")
        _run(f"chmod +x {tool_path}")
    return f"./{tool_path}"

def fetch_ccres(chrom, start, end, out_bed="data/cache/ccre_region.bed"):
    os.makedirs(os.path.dirname(out_bed), exist_ok=True)
    tool = fetch_bigbedtobed_tool()
    remote_bb = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb"
    cmd = f"{tool} {remote_bb} -chrom={to_ucsc_chrom(chrom)} -start={start} -end={end} {out_bed}"
    r = _run(cmd, check=False)
    if r.returncode != 0 or not os.path.exists(out_bed) or os.path.getsize(out_bed) == 0:
        print("WARNING: cCRE fetch failed or returned empty for this region")
        return pd.DataFrame(columns=["Chromosome", "Start", "End", "ccre_label"])
    
    cols = ["Chromosome", "Start", "End", "name", "score", "strand", "thickStart", "thickEnd",
            "reserved", "ccre_group", "ccre_group2", "zscore", "ucscLabel", "accession", "ccre_label"]
    df = pd.read_csv(out_bed, sep="\\t", header=None, names=cols)
    df["Chromosome"] = to_plain_chrom(chrom)
    return df[["Chromosome", "Start", "End", "ccre_label"]]

def flag_interval_overlap(variants_df, reg_intervals_df):
    variants_df = variants_df.reset_index(drop=True).copy()
    variants_df["_row_id"] = variants_df.index

    if reg_intervals_df.empty:
        variants_df["interval_evidence"] = [[] for _ in range(len(variants_df))]
        return variants_df.drop(columns=["_row_id"])

    var_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": variants_df["chrom"],
        "Start": variants_df["pos"] - 1,
        "End": variants_df["pos"],
        "row_id": variants_df["_row_id"],
    }))
    reg_pr = pr.PyRanges(reg_intervals_df)

    joined = var_pr.join(reg_pr).df
    if joined.empty or "row_id" not in joined.columns:
        variants_df["interval_evidence"] = [[] for _ in range(len(variants_df))]
        return variants_df.drop(columns=["_row_id"])

    evidence_map = joined.groupby("row_id")["ccre_label"].apply(list).to_dict()
    variants_df["interval_evidence"] = variants_df["_row_id"].map(lambda i: evidence_map.get(i, []))
    return variants_df.drop(columns=["_row_id"])
""",

    "src/regvar/layer0/vep.py": """import time
import requests
import pandas as pd

def vep_annotate_batch(variants_df, batch_size=200, timeout=30, max_retries=3):
    results = []
    variants_list = variants_df[["chrom", "pos", "ref", "alt"]].to_dict("records")
    n_requested = len(variants_list)

    for i in range(0, n_requested, batch_size):
        batch = variants_list[i:i + batch_size]
        variant_strings = [f"{v['chrom']} {v['pos']} . {v['ref']} {v['alt']} . . ." for v in batch]

        attempt = 0
        while attempt <= max_retries:
            try:
                r = requests.post(
                    "https://rest.ensembl.org/vep/human/region",
                    headers={"Content-Type": "application/json", "Accept": "application/json"},
                    json={"variants": variant_strings, "regulatory": 1},
                    timeout=timeout,
                )
                if r.status_code == 200:
                    results.extend(r.json())
                    break
                elif r.status_code in (429, 503):
                    wait = 2 ** attempt
                    print(f"Batch {i}: rate-limited ({r.status_code}), retrying in {wait}s")
                    time.sleep(wait)
                    attempt += 1
                else:
                    print(f"Batch {i} failed permanently: {r.status_code}, {r.text[:200]}")
                    break
            except requests.exceptions.RequestException as e:
                wait = 2 ** attempt
                print(f"Batch {i}: {type(e).__name__}, retrying in {wait}s")
                time.sleep(wait)
                attempt += 1
        time.sleep(1)

    return results

def parse_vep_results(vep_results, canonical_transcript, target_gene_id):
    rows = []
    for res in vep_results:
        pos = res.get("start")
        chrom = res.get("seq_region_name")
        allele = res.get("allele_string")
        if not allele or "/" not in allele:
            continue
        ref, alt = allele.split("/")[:2]

        gene_restricted_terms = set()
        canonical_terms = []
        for tc in res.get("transcript_consequences", []):
            if tc.get("gene_id") == target_gene_id:
                gene_restricted_terms.update(tc.get("consequence_terms", []))
                if tc.get("transcript_id") == canonical_transcript:
                    canonical_terms = tc.get("consequence_terms", [])

        reg_terms = set()
        for rfc in res.get("regulatory_feature_consequences", []):
            reg_terms.update(rfc.get("consequence_terms", []))

        rows.append({
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "vep_canonical_terms": canonical_terms,
            "vep_all_transcript_terms": sorted(gene_restricted_terms),
            "vep_regulatory_terms": sorted(reg_terms),
        })

    return pd.DataFrame(rows)
""",

    "src/regvar/layer0/pipeline.py": """import os
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
""",

    "src/regvar/cli.py": """import argparse
from regvar.layer0.pipeline import run_layer0

def main():
    parser = argparse.ArgumentParser(description="Regulatory Variant Prioritization Framework")
    parser.add_argument("command", choices=["run"], help="Command to execute")
    parser.add_argument("--layer", type=int, required=True, help="Layer to run (e.g., 0)")
    parser.add_argument("--config", type=str, required=True, help="Path to config YAML")

    args = parser.parse_args()

    if args.command == "run":
        if args.layer == 0:
            print(f"Running Layer 0 with config {args.config}...")
            run_layer0(args.config)
        else:
            print(f"Layer {args.layer} is not yet implemented.")

if __name__ == "__main__":
    main()
"""
}

for path, content in files.items():
    full_path = os.path.join(base, path)
    os.makedirs(os.path.dirname(full_path), exist_ok=True)
    with open(full_path, "w", encoding="utf-8") as f:
        f.write(content)

print("Layer 0 python files generated.")
