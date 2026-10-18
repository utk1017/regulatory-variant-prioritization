# Framework Overview

The Regulatory Variant Prioritization Framework operates in a multi-layered architecture:

- **Layer 0 (Variant Intake & Harmonization):** Flags whether a variant sits in a regulatory region (cCREs, promoters, UTRs) using interval overlap and Ensembl VEP annotations. Outputs the candidate truth-set.
- **Layer 1 (Multi-Source Evidence Collection):** Extracts and standardizes observations across gnomAD (population AF), GWAS Catalog (associations), GTEx Portal (eQTLs), MPRA, and CRISPR screens.
- **Layer 2 (Evidence Rules & Completability):** Determines which `absent` modalities trigger computational completion models based on modality metadata (`completable=True`).
- **Layer 3 (AlphaGenome):** Computational completion for missing experimental modalities. *(In progress)*
- **Layer 5 (Bayesian Evidence Fusion):** Combines multi-source experimental and computational evidence into a calibrated regulatory posterior score. *(In progress)*
- **Layer 6 (Calibration & Stratification):** Calibrates prior and likelihood distributions against ClinVar benchmarks across evidence strata. *(In progress)*
- **Layer 7 (Variant Prioritization & Output):** Produces the final ranked variant predictions with full provenance manifests. *(In progress)*

See `README.md` for Layer 0 and Layer 1 execution details.
