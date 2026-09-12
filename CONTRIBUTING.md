# Contributing to the Regulatory Variant Framework

Welcome! To maintain the integrity of this scientific pipeline, please adhere to the following strict design rules when contributing.

## 1. `schema.py` is the load-bearing file
Every layer produces `EvidenceRecord`-shaped output. This schema is defined once in `src/regvar/common/schema.py` and imported everywhere. Do not invent new column names when adding layers. 

## 2. Provenance is data, not a side effect
Every pipeline run writes a `manifest.json`. If a result can't be traced back to the exact inputs (snapshot date, genome build, config file) that produced it, treat that as a bug.

## 3. Genes are config, not code
Adding a new gene or expanding to a panel should never require touching the `src/` directory. All gene definitions live in `configs/`.

## 4. `data/` is fully gitignored
Never commit reference FASTAs, ClinVar snapshots, or generated CSVs to the repository. If you need the actual data versioned, use DVC or git-lfs — don't fight git with raw CSVs.

## 5. Notebooks stay thin
Logic lives in `src/regvar/`, tested independently. Notebooks under `notebooks/` are demo wrappers that call the installed package. 

## 6. Layer 0's scope boundary is deliberate
Layer 0 answers "is this variant a regulatory candidate, and what's the evidence". It does NOT answer "which kind of element is this" or "does it break the mechanism". Scope creep will be rejected.
