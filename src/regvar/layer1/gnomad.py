"""gnomAD population allele frequency collection for Layer 1.

Annotation-only modality: gnomAD is completable=False (AlphaGenome cannot
predict allele frequency) and use_in_fusion=False (keeping AF out of the
Bayesian posterior prevents label-calibration leakage in Layer 6, where the
benign set is defined by allele frequency).

The MAF derived here is also used by the eQTL module to decide whether a
non-hit is a real null (common variant, testable in GTEx) or merely absent
(rare variant, never tested by GTEx).

Implementation notes
--------------------
- Region queries are windowed from the candidate positions using
  :func:`~regvar.layer1.variants.plan_windows` so distal candidates such as
  LCR variants are covered by construction.
- Only ``ac`` and ``an`` are requested; AF is computed locally (``ac / an``).
- A failed window marks the variants inside it ``absent / query_failed``
  rather than silently treating them as unseen.
- ``gnomad_request_delay_s`` only sleeps *between real network calls*; cache
  hits return instantly without consuming any inter-request budget.
"""
from __future__ import annotations

import time
from collections import defaultdict
from typing import Any, Callable, Dict, List, Optional

from regvar.layer1.schema import make_record
from regvar.layer1.variants import Variant, canon_key, plan_windows


GNOMAD_API = "https://gnomad.broadinstitute.org/api"

GNOMAD_REGION_QUERY = """
query($chrom: String!, $start: Int!, $stop: Int!, $dataset: DatasetId!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    variants(dataset: $dataset) {
      variant_id
      pos
      ref
      alt
      genome { ac an }
    }
  }
}
"""


def collect_gnomad_evidence(
    variants: List[Variant],
    dataset: str,
    window_bp: int,
    delay_s: float,
    fetch: Callable,
) -> Dict[str, Any]:
    """Query gnomAD for all candidate variants and return a result dict.

    The result dict has the form::

        {
            "windows": [{"chrom", "start", "stop", "error", "n_variants", "from_cache"}, ...],
            "by_key":  {canon_key: {"ac", "an", "af", "variant_id"}, ...},
        }
    """
    result: Dict[str, Any] = {"windows": [], "by_key": {}}
    by_chrom: dict = defaultdict(list)
    for v in variants:
        by_chrom[v.chrom].append(v.pos)

    first_network_call = True
    for chrom, positions in by_chrom.items():
        for s, e in plan_windows(positions, window_bp):
            start, stop = max(1, s - 1), e + 1
            payload, err, from_cache = fetch(
                "POST",
                GNOMAD_API,
                json_body={
                    "query": GNOMAD_REGION_QUERY,
                    "variables": {
                        "chrom": chrom,
                        "start": start,
                        "stop": stop,
                        "dataset": dataset,
                    },
                },
            )
            # Apply rate-limit delay only for real network calls, not cache hits.
            if not from_cache and delay_s:
                if not first_network_call:
                    time.sleep(delay_s)
                first_network_call = False

            if not err and payload and payload.get("errors"):
                err = "graphql_error: " + str(payload["errors"][0].get("message", ""))[:120]

            win: Dict[str, Any] = {
                "chrom": chrom,
                "start": s,
                "stop": e,
                "error": err,
                "n_variants": 0,
                "from_cache": from_cache,
            }
            if not err:
                region = ((payload or {}).get("data") or {}).get("region") or {}
                for gv in region.get("variants") or []:
                    genome = gv.get("genome")
                    if not genome or not genome.get("an") or not genome.get("ac"):
                        continue
                    key = canon_key(chrom, gv["pos"], gv["ref"], gv["alt"])
                    result["by_key"][key] = {
                        "ac": genome["ac"],
                        "an": genome["an"],
                        "af": genome["ac"] / genome["an"],
                        "variant_id": gv.get("variant_id"),
                    }
                    win["n_variants"] += 1
            result["windows"].append(win)
    return result


def build_gnomad_records(gnomad: Dict[str, Any], variants: List[Variant], dataset: str):
    """Convert the *gnomad* result dict into a list of :class:`~regvar.layer1.schema.EvidenceRecord`."""
    records = []
    for v in variants:
        hit = gnomad["by_key"].get(v.key)
        if hit:
            records.append(
                make_record(
                    v,
                    "gnomad_af",
                    "curated_clinical",
                    "observed_positive",
                    context=f"{dataset}:genomes",
                    af=hit["af"],
                    ac=hit["ac"],
                    an=hit["an"],
                    gnomad_variant_id=hit["variant_id"],
                )
            )
            continue
        failed = [
            w
            for w in gnomad["windows"]
            if w["chrom"] == v.chrom
            and w["start"] <= v.pos <= w["stop"]
            and w["error"]
        ]
        if failed:
            records.append(
                make_record(
                    v,
                    "gnomad_af",
                    "curated_clinical",
                    "absent",
                    context=f"{dataset}:genomes",
                    reason="query_failed",
                    error=failed[0]["error"],
                )
            )
        else:
            records.append(
                make_record(
                    v,
                    "gnomad_af",
                    "curated_clinical",
                    "absent",
                    context=f"{dataset}:genomes",
                    reason="not_observed_in_genome_callset",
                )
            )
    return records


def maf_from_gnomad(gnomad: Dict[str, Any]) -> Dict:
    """Derive a ``{canon_key: maf}`` lookup from gnomAD results.

    MAF = min(AF, 1-AF) so that heterozygous alternate and reference alleles
    are treated symmetrically.  Used by the eQTL module.
    """
    return {
        key: min(hit["af"], 1 - hit["af"])
        for key, hit in gnomad["by_key"].items()
    }
