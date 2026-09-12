import os
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
