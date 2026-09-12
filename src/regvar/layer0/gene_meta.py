import requests

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
