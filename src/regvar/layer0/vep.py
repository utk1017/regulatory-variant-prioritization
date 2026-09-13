import time
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

    n_returned = len(results)
    if n_returned < n_requested:
        print(f"WARNING: requested {n_requested} variants, got {n_returned} VEP results "
              f"({n_requested - n_returned} lost to failed/skipped batches)")
    return results

def parse_vep_results(vep_results, canonical_transcript, target_gene_id):
    rows = []
    n_missing_allele = 0
    for res in vep_results:
        pos = res.get("start")
        chrom = res.get("seq_region_name")
        allele = res.get("allele_string")
        if not allele or "/" not in allele:
            n_missing_allele += 1
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

    if n_missing_allele:
        print(f"NOTE: {n_missing_allele} VEP results had no parseable allele_string and were "
              f"dropped (these will not match any variant downstream)")
    return pd.DataFrame(rows)
