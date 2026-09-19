# Framework Overview

The Regulatory Variant Prioritization Framework operates in a multi-layered architecture:

- **Layer 0 (Variant Intake):** Flags whether a variant sits in a regulatory region (cCREs, promoters, UTRs) using interval overlap and Ensembl VEP annotations. Outputs the candidate truth-set.
- **Layer 1 (Experimental Evidence):** Incorporates GWAS, eQTL, MPRA, and CRISPR screens. *(In progress)*
- **Layer 3 (AlphaGenome):** Computational completion for missing modalities. *(In progress)*
- **Layer 7 (Integration):** Final Bayesian integration against ClinVar reference sets. *(In progress)*

See `README.md` for Layer 0 execution details.
