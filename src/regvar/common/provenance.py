"""Implements CONTRIBUTING.md rule #2: provenance is data, not a side effect.

Every pipeline run writes a manifest.json capturing exactly what produced its output --
snapshot dates, genome build, config used, git commit -- so a result printed to stdout and then
lost in a scrollback isn't the only record of how it was generated.
"""
import json
import os
import subprocess
import hashlib
from datetime import datetime, timezone


def _git_commit_hash():
    try:
        r = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True, timeout=5
        )
        return r.stdout.strip()
    except Exception:
        return None


def _file_hash(path, algo="sha256"):
    if not path or not os.path.exists(path):
        return None
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(
    out_dir,
    *,
    gene_symbol,
    gene_meta,
    clinvar_snapshot_date,
    clinvar_date_was_pinned,
    ref_fasta_path,
    config_path,
    n_variants_total,
    n_variants_regulatory,
    extra=None,
):
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit_hash(),
        "gene_symbol": gene_symbol,
        "gene_metadata": gene_meta,
        "clinvar_snapshot_date": clinvar_snapshot_date,
        "clinvar_date_was_pinned": clinvar_date_was_pinned,
        "reference_genome_build": "GRCh38/hg38",
        "reference_fasta_path": ref_fasta_path,
        "config_path": config_path,
        "config_sha256": _file_hash(config_path),
        "n_variants_total": n_variants_total,
        "n_variants_regulatory": n_variants_regulatory,
    }
    if not clinvar_date_was_pinned:
        manifest["warning"] = (
            "clinvar_snapshot_date was auto-discovered, not pinned in the config. "
            "Copy this date into the config's clinvar_snapshot_date field for reproducible re-runs."
        )
    if extra:
        manifest.update(extra)

    os.makedirs(out_dir, exist_ok=True)
    manifest_path = os.path.join(out_dir, "manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    return manifest_path
