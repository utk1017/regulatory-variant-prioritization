import os
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_ALIGN_VERTICAL
from docx.oxml import parse_xml, OxmlElement
from docx.oxml.ns import nsdecls, qn

def set_cell_background(cell, fill_hex):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = parse_xml(f'<w:shd {nsdecls("w")} w:fill="{fill_hex}"/>')
    tcPr.append(shd)

def set_cell_padding(cell, top=120, bottom=120, left=180, right=180):
    tcPr = cell._tc.get_or_add_tcPr()
    tcMar = parse_xml(
        f'<w:tcMar {nsdecls("w")}>'
        f'<w:top w:w="{top}" w:type="dxa"/>'
        f'<w:bottom w:w="{bottom}" w:type="dxa"/>'
        f'<w:left w:w="{left}" w:type="dxa"/>'
        f'<w:right w:w="{right}" w:type="dxa"/>'
        f'</w:tcMar>'
    )
    tcPr.append(tcMar)

def style_table(table, col_widths=None, header_bg="1E3A8A", alt_bg="F8FAFC"):
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, row in enumerate(table.rows):
        # CantSplit rule to prevent row breaks across pages
        trPr = row._tr.get_or_add_trPr()
        trPr.append(parse_xml(f'<w:cantSplit {nsdecls("w")}/>'))
        
        # Header repeat on every page
        if i == 0:
            trPr.append(parse_xml(f'<w:tblHeader {nsdecls("w")}/>'))
            
        for j, cell in enumerate(row.cells):
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
            set_cell_padding(cell, top=140, bottom=140, left=180, right=180)
            if i == 0:
                set_cell_background(cell, header_bg)
                for p in cell.paragraphs:
                    p.paragraph_format.space_before = Pt(2)
                    p.paragraph_format.space_after = Pt(2)
                    for run in p.runs:
                        run.font.bold = True
                        run.font.color.rgb = RGBColor(255, 255, 255)
                        run.font.name = "Calibri"
                        run.font.size = Pt(10)
            else:
                if i % 2 == 1:
                    set_cell_background(cell, "FFFFFF")
                else:
                    set_cell_background(cell, alt_bg)
                for p in cell.paragraphs:
                    p.paragraph_format.space_before = Pt(2)
                    p.paragraph_format.space_after = Pt(2)
                    for run in p.runs:
                        run.font.name = "Calibri"
                        run.font.size = Pt(9.5)
                        run.font.color.rgb = RGBColor(30, 41, 59)
            if col_widths and j < len(col_widths):
                cell.width = col_widths[j]

def add_callout(doc, text, title="KEY ARCHITECTURAL DECISION"):
    tbl = doc.add_table(rows=1, cols=1)
    tbl.alignment = WD_TABLE_ALIGNMENT.CENTER
    cell = tbl.cell(0, 0)
    cell.width = Inches(6.5)
    set_cell_background(cell, "F0F4F8")
    set_cell_padding(cell, top=160, bottom=160, left=240, right=200)
    
    # Left border only (thick navy border)
    tcPr = cell._tc.get_or_add_tcPr()
    borders = parse_xml(
        f'<w:tcBorders {nsdecls("w")}>'
        f'<w:top w:val="none"/>'
        f'<w:left w:val="single" w:sz="36" w:space="0" w:color="1E3A8A"/>'
        f'<w:bottom w:val="none"/>'
        f'<w:right w:val="none"/>'
        f'</w:tcBorders>'
    )
    tcPr.append(borders)
    
    p = cell.paragraphs[0]
    p.paragraph_format.space_before = Pt(2)
    p.paragraph_format.space_after = Pt(4)
    run_t = p.add_run(f"[{title}] ")
    run_t.font.bold = True
    run_t.font.size = Pt(10)
    run_t.font.color.rgb = RGBColor(30, 58, 138)
    
    run_b = p.add_run(text)
    run_b.font.size = Pt(9.5)
    run_b.font.color.rgb = RGBColor(51, 65, 85)
    run_b.font.name = "Calibri"
    
    # Add an empty paragraph after callout for breathing room
    p_after = doc.add_paragraph()
    p_after.paragraph_format.space_before = Pt(0)
    p_after.paragraph_format.space_after = Pt(6)

def build_word_document():
    doc = Document()
    
    # Page setup - 1 inch margins
    sections = doc.sections
    for section in sections:
        section.top_margin = Inches(1.0)
        section.bottom_margin = Inches(1.0)
        section.left_margin = Inches(1.0)
        section.right_margin = Inches(1.0)

    # -------------------------------------------------------------
    # Cover / Header Title
    # -------------------------------------------------------------
    p_title = doc.add_paragraph()
    p_title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p_title.paragraph_format.space_before = Pt(12)
    p_title.paragraph_format.space_after = Pt(4)
    run_main_title = p_title.add_run("Regulatory Variant Prioritization Framework")
    run_main_title.font.name = "Calibri"
    run_main_title.font.size = Pt(24)
    run_main_title.font.bold = True
    run_main_title.font.color.rgb = RGBColor(30, 58, 138) # Deep Navy

    p_sub = doc.add_paragraph()
    p_sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p_sub.paragraph_format.space_before = Pt(0)
    p_sub.paragraph_format.space_after = Pt(14)
    run_sub = p_sub.add_run("Technical Architecture & Progress Report: Layers 0 to 2")
    run_sub.font.name = "Calibri"
    run_sub.font.size = Pt(15)
    run_sub.font.color.rgb = RGBColor(71, 85, 105) # Slate

    # Metadata Card
    meta_tbl = doc.add_table(rows=2, cols=2)
    meta_tbl.alignment = WD_TABLE_ALIGNMENT.CENTER
    style_table(meta_tbl, col_widths=[Inches(3.25), Inches(3.25)], header_bg="0F172A", alt_bg="F1F5F9")
    meta_tbl.rows[0].cells[0].paragraphs[0].text = "Project Team & Investigators"
    meta_tbl.rows[0].cells[1].paragraphs[0].text = "Project & Milestones"
    
    p0 = meta_tbl.rows[1].cells[0].paragraphs[0]
    p0.text = "Utkarsh Gupta & Kamakshi Mudgal\nAdvanced Genomic Computing / Bio-AI"
    p1 = meta_tbl.rows[1].cells[1].paragraphs[0]
    p1.text = "Repository: regvar | Build Date: November 2026\nStatus: Layers 0–2 Verified | Layers 3–7 In Progress"
    style_table(meta_tbl, col_widths=[Inches(3.25), Inches(3.25)], header_bg="1E3A8A", alt_bg="F8FAFC")

    p_sep = doc.add_paragraph()
    p_sep.paragraph_format.space_before = Pt(10)
    p_sep.paragraph_format.space_after = Pt(10)

    # -------------------------------------------------------------
    # 1. Executive Summary & Problem Formulation
    # -------------------------------------------------------------
    h1 = doc.add_heading("1. Executive Summary & Problem Formulation", level=1)
    h1.paragraph_format.space_before = Pt(14)
    h1.paragraph_format.space_after = Pt(6)

    doc.add_paragraph(
        "Genome-wide association studies (GWAS) and whole-genome sequencing (WGS) have consistently demonstrated "
        "that over 90% of disease-associated genetic variants lie within non-coding regulatory sequences. "
        "However, clinical interpretation frameworks (such as ACMG/AMP standards) remain overwhelmingly optimized "
        "for protein-coding regions, where amino acid substitutions, frameshifts, and nonsense mutations produce direct "
        "biochemical consequences. Interpreting non-coding single-nucleotide variants (SNVs) requires synthesizing noisy, "
        "sparse, and multimodal regulatory evidence across population genetics, epigenomics, and deep sequence foundation models."
    )

    doc.add_paragraph(
        "Current clinical and bioinformatic pipelines suffer from four critical failure modes:\n"
        "• Severe Experimental Sparsity: High-throughput functional assays (MPRA, CRISPRi, eQTL) experimentally test fewer than 1% of possible non-coding variants.\n"
        "• Confounding Missingness vs. Negative Measurements: Existing algorithms frequently treat an unmeasured variant (a coverage gap) as equivalent to a variant tested and proven null (an experimental negative). This conflation severely distorts Bayesian evidence weighting.\n"
        "• Benchmark Contamination & Circularity: Evaluated computational models are frequently trained on the same variant sets used for clinical calibration, causing severe optimistic bias.\n"
        "• Locus & Gene Cluster Contamination: In dense genomic clusters (such as the human beta-globin locus on chromosome 11), naïve spatial proximity misattributes regulatory effects to neighboring genes and pseudogenes."
    )

    add_callout(
        doc,
        "Our framework resolves these fundamental challenges by combining: (1) curated clinical ground truth (ClinVar >= 1 star), "
        "(2) multi-source experimental evidence across five modalities, (3) a 3-state missingness taxonomy disambiguating biological gaps from experimental nulls, "
        "(4) zero-shot foundation model inference via AlphaGenome, and (5) evidence-stratified Bayesian calibration.",
        title="FRAMEWORK MISSION"
    )

    # -------------------------------------------------------------
    # 2. System Architecture
    # -------------------------------------------------------------
    h2 = doc.add_heading("2. End-to-End System Architecture (Layers 0 to 7)", level=1)
    h2.paragraph_format.space_before = Pt(14)
    h2.paragraph_format.space_after = Pt(6)

    doc.add_paragraph(
        "The Regulatory Variant Prioritization Framework (regvar) is organized into eight decoupled, modular layers. "
        "Each layer communicates via strictly validated schema contracts and emits an immutable manifest containing "
        "full provenance metadata, cryptographic file hashes, and git commit IDs."
    )

    layers_tbl = doc.add_table(rows=9, cols=4)
    layers_tbl.rows[0].cells[0].paragraphs[0].text = "Layer"
    layers_tbl.rows[0].cells[1].paragraphs[0].text = "Name & Scope"
    layers_tbl.rows[0].cells[2].paragraphs[0].text = "Core Responsibilities"
    layers_tbl.rows[0].cells[3].paragraphs[0].text = "Status"

    layers_data = [
        ("Layer 0", "Variant Intake & Harmonization", "Coordinates, ClinVar snapshot pinning, ENCODE cCRE overlap, gene-restricted VEP union", "Completed & Tested"),
        ("Layer 1", "Multi-Source Evidence Collection", "Standardizes gnomAD v4, GWAS Catalog, GTEx v8, MPRA, and CRISPR into 3-state missingness matrix", "Completed & Tested"),
        ("Layer 2", "Gap Detection & Trigger Rules", "Disambiguates biological gaps vs tested nulls; compiles AlphaGenome trigger queue & agreement benchmark", "Completed & Tested"),
        ("Layer 3", "AlphaGenome Completion", "Foundation model execution across 1 Mb sequence contexts for missing completable modalities", "In Progress"),
        ("Layer 4", "Computational Feature Extraction", "Computes allelic delta tracks (expression, chromatin, TF binding) and harmonizes proxy quality", "Planned"),
        ("Layer 5", "Bayesian Evidence Fusion", "Hierarchical multimodal likelihood integration with proxy downweighting", "Planned"),
        ("Layer 6", "Calibration & Stratification", "ACMG-style Bayesian calibration against ClinVar reference sets across Tiers 1-3", "Planned"),
        ("Layer 7", "Clinical Output & Prioritization", "Ranked variant reports, mechanistic audit manifests, and clinical classifications", "Planned"),
    ]

    for idx, (l_id, l_name, l_resp, l_stat) in enumerate(layers_data, start=1):
        layers_tbl.rows[idx].cells[0].paragraphs[0].text = l_id
        layers_tbl.rows[idx].cells[1].paragraphs[0].text = l_name
        layers_tbl.rows[idx].cells[2].paragraphs[0].text = l_resp
        layers_tbl.rows[idx].cells[3].paragraphs[0].text = l_stat

    style_table(layers_tbl, col_widths=[Inches(0.9), Inches(1.8), Inches(2.6), Inches(1.2)], header_bg="1E3A8A")

    # -------------------------------------------------------------
    # 3. Layer 0 Deep Dive
    # -------------------------------------------------------------
    h3 = doc.add_heading("3. Deep Dive: Layer 0 — Variant Intake & Harmonization", level=1)
    h3.paragraph_format.space_before = Pt(16)
    h3.paragraph_format.space_after = Pt(6)

    doc.add_paragraph(
        "Scope Boundary: Layer 0 answers exactly one question: 'Does a candidate variant sit within a regulatory region of the target gene, "
        "and what exact genomic evidence supports that call?' It does not attempt to predict element mechanism or calculate scores."
    )

    doc.add_paragraph(
        "Key Methodological Components:\n"
        "1. Dynamic Coordinate & Strand-Aware Promoter Windows: Automatically resolves chromosome, coordinates, and strand from Ensembl. "
        "Computes a strand-aware promoter interval (default: [-2500, 0] bp relative to TSS) and a flanking boundary (default: +/- 100 kb) "
        "to capture distal elements such as the Locus Control Region (LCR).\n"
        "2. Pinned ClinVar Snapshotting: To guarantee scientific reproducibility, the pipeline auto-discovers the latest snapshot date once, "
        "and pins it into the gene YAML config (e.g. clinvar_snapshot_date: '20240915'). Variants are extracted via bcftools and normalized "
        "against GRCh38 using bcftools norm -m -both and deduplicated.\n"
        "3. Review-Star Quality Filtering: Excludes conflicting classifications and unasserted entries, enforcing review_stars >= 1 "
        "(single submitter with criteria, expert panel, or practice guideline).\n"
        "4. Interval-First Regulatory Classification: Downloads ENCODE cCREs (v4) via UCSC bigBedToBed and performs vectorized interval overlap "
        "joins via PyRanges. Captures proximal enhancer-like signatures (pELS), distal enhancer-like signatures (dELS), and computed promoters.\n"
        "5. Gene-Restricted Transcript Consequence Union: Queries Ensembl VEP REST API in batches of 200 variants. Uniquely restricts transcript consequences "
        "to the target gene's Ensembl Gene ID while unioning terms across all isoforms."
    )

    l0_edge_tbl = doc.add_table(rows=4, cols=3)
    l0_edge_tbl.rows[0].cells[0].paragraphs[0].text = "Edge Case / Failure Mode"
    l0_edge_tbl.rows[0].cells[1].paragraphs[0].text = "Biological & Algorithmic Cause"
    l0_edge_tbl.rows[0].cells[2].paragraphs[0].text = "Solution & Regression Lock"

    l0_edge_data = [
        ("Distal Enhancer (LCR) Dropping", "Variants in distal enhancers (e.g., HBB LCR at chr11:5,280,000, 50 kb upstream) have zero transcript overlap and are labeled intergenic by VEP.", "Interval-First Rule: Overlap with an ENCODE cCRE flags a variant as regulatory regardless of VEP transcript terms. Verified by test_lcr_survival.py."),
        ("Neighboring Gene Contamination", "In dense gene clusters (e.g. HBB alongside HBD and HBBP1 pseudogene), wide-window VEP queries return consequences belonging to neighbor genes.", "Gene-ID Restriction: Transcript consequences are strictly filtered to the target Ensembl Gene ID (ENSG00000244734). Verified by test_gene_scoped_vep.py."),
        ("Alternative TSS Disagreement", "Pathogenic 5' UTR variant (chr11:5,227,172 G>A) was misclassified as upstream because the canonical transcript's TSS was offset from an alternative transcript's TSS.", "All-Transcript Union: Consequence terms are unioned across all transcripts of the gene, correctly capturing 5_prime_UTR_variant. Verified by test_canonical_union.py."),
    ]
    for idx, (ec, cause, sol) in enumerate(l0_edge_data, start=1):
        l0_edge_tbl.rows[idx].cells[0].paragraphs[0].text = ec
        l0_edge_tbl.rows[idx].cells[1].paragraphs[0].text = cause
        l0_edge_tbl.rows[idx].cells[2].paragraphs[0].text = sol

    style_table(l0_edge_tbl, col_widths=[Inches(1.8), Inches(2.2), Inches(2.5)], header_bg="1E3A8A")

    # -------------------------------------------------------------
    # 4. Layer 1 Deep Dive
    # -------------------------------------------------------------
    h4 = doc.add_heading("4. Deep Dive: Layer 1 — Multi-Source Evidence Collection", level=1)
    h4.paragraph_format.space_before = Pt(16)
    h4.paragraph_format.space_after = Pt(6)

    doc.add_paragraph(
        "Layer 1 ingests the candidate truth-set from Layer 0 and collects empirical evidence across five distinct biological modalities. "
        "Crucially, it replaces ambiguous null values with an explicit 3-state missingness taxonomy:"
    )

    doc.add_paragraph(
        "• observed_positive: Statistically significant experimental signal (GWAS p <= 5e-8, eQTL p <= 0.05, significant MPRA skew).\n"
        "• observed_null: The variant was testable and assayed, but showed no statistically significant effect (e.g., common variant MAF >= 0.01 in GTEx Whole Blood with no expression alteration).\n"
        "• absent: The variant was not testable, unobserved in the catalog, or affected by a query failure. Always accompanied by a structured reason code."
    )

    add_callout(
        doc,
        "Crucial Biological Distinction: A rare beta-thalassemia pathogenic mutation unmeasured in GTEx because of statistical power (MAF < 0.01) "
        "is marked 'absent' (reason: not_testable_rare_or_unknown_frequency). It must NEVER be marked 'observed_null', which would wrongly penalize it "
        "as an experimental negative.",
        title="CLINICAL SENSITIVITY PRINCIPLE"
    )

    doc.add_heading("Modality Governance & Policies (MODALITY_META)", level=2)
    mod_tbl = doc.add_table(rows=6, cols=5)
    mod_tbl.rows[0].cells[0].paragraphs[0].text = "Modality"
    mod_tbl.rows[0].cells[1].paragraphs[0].text = "Data Source & Scope"
    mod_tbl.rows[0].cells[2].paragraphs[0].text = "Completable by AlphaGenome?"
    mod_tbl.rows[0].cells[3].paragraphs[0].text = "Use in Fusion?"
    mod_tbl.rows[0].cells[4].paragraphs[0].text = "Primary Function"

    mod_data = [
        ("gwas", "GWAS Catalog REST API v2", "False (Population phenotype)", "True", "Locus/trait association"),
        ("eqtl", "GTEx Portal REST API v2", "True (Direct expression track)", "True", "Tissue expression alteration"),
        ("mpra", "Curated Local Tables", "True (Chromatin / initiation proxy)", "True", "Episomal reporter skew"),
        ("crispr", "Curated Local Tables", "True (Expression perturbation proxy)", "True", "Endogenous perturbation"),
        ("gnomad_af", "gnomAD v4 GraphQL API", "False (Cohort property)", "False (Stratification only)", "MAF proxy / benign leaning"),
    ]
    for idx, (m, src, comp, fus, func) in enumerate(mod_data, start=1):
        mod_tbl.rows[idx].cells[0].paragraphs[0].text = m
        mod_tbl.rows[idx].cells[1].paragraphs[0].text = src
        mod_tbl.rows[idx].cells[2].paragraphs[0].text = comp
        mod_tbl.rows[idx].cells[3].paragraphs[0].text = fus
        mod_tbl.rows[idx].cells[4].paragraphs[0].text = func

    style_table(mod_tbl, col_widths=[Inches(1.0), Inches(1.8), Inches(1.4), Inches(1.0), Inches(1.3)], header_bg="1E3A8A")

    doc.add_paragraph(
        "Robust Engineering in Layer 1:\n"
        "• Resilient HTTP Cache: Disk-backed cache with SHA-256 keys and 7-day TTL. Rejects HTTP 200 payloads carrying GraphQL errors.\n"
        "• Dynamic Rate-Limiting: 6-second sleep enforced on live gnomAD calls, bypassed instantly on cache hits.\n"
        "• Windowed GWAS Error Scoping: Variant positions clustered into +/- 25 kb windows; errors in one window do not contaminate adjacent windows.\n"
        "• Complete Matrix Enforcement: validate_matrix_complete asserts every variant has a record for every modality."
    )

    # -------------------------------------------------------------
    # 5. Layer 2 Deep Dive
    # -------------------------------------------------------------
    h5 = doc.add_heading("5. Deep Dive: Layer 2 — Evidence Rules, Missingness & Trigger Compilation", level=1)
    h5.paragraph_format.space_before = Pt(16)
    h5.paragraph_format.space_after = Pt(6)

    doc.add_paragraph(
        "Layer 2 is the routing and rule engine bridging Layer 1 empirical evidence with Layer 3 AlphaGenome execution and Layer 6 calibration. "
        "It performs five central operations:"
    )

    doc.add_paragraph(
        "1. Semantic Reason Disambiguation Engine: Maps evidence observations into five operationally distinct categories:\n"
        "   - POSITIVE: Real biological hit -> routed to Agreement Benchmark Queue.\n"
        "   - NULL_EXPERIMENTAL: Real experimental negative -> routed to Agreement Benchmark Queue. NEVER completed by AlphaGenome.\n"
        "   - GAP_BIOLOGICAL: True coverage gap -> eligible for AlphaGenome trigger queue.\n"
        "   - ERROR_TRANSIENT: API or network failure -> flagged as warning, blocked from AlphaGenome completion.\n"
        "   - UNKNOWN_REASON: Unrecognized reason string -> flagged as warning, blocked from AlphaGenome completion.\n"
        "2. Completability Filtering: AlphaGenome foundation models can predict functional tracks (gene expression, chromatin accessibility, transcription initiation), "
        "but cannot predict population statistics (GWAS association or gnomAD allele frequency). Non-completable modalities are strictly prohibited from entering the trigger queue.\n"
        "3. AlphaGenome Trigger Queue Compilation: Emits layer2_alphagenome_triggers.csv with target gene, target tissue/cell context, "
        "receptive field (1 Mb for expression, 2 kb for local reporter), and output track mapping.\n"
        "4. Segregated Agreement Benchmark Queue: Real experimental measurements (observed_positive and observed_null) on completable modalities are isolated "
        "into layer2_agreement_checks.csv (validation_purpose='internal_agreement_benchmark_only'). Used solely to benchmark AlphaGenome zero-shot accuracy "
        "without circular data leakage into Bayesian fusion.\n"
        "5. Evidence Stratification (Tiers 1 to 3): Computes the count of tested allele-resolved experimental modalities N_exp in {eqtl, mpra, crispr}:\n"
        "   - Tier 1 (Evidence-Rich): N_exp >= 2 modalities tested.\n"
        "   - Tier 2 (Moderate Evidence): N_exp = 1 modality tested.\n"
        "   - Tier 3 (Evidence-Sparse / Zero-Coverage): N_exp = 0 modalities tested. (Primary beneficiary of AlphaGenome completion)."
    )

    add_callout(
        doc,
        "Safety Guarantee: Why Experimental Nulls Never Trigger AlphaGenome: If an MPRA assay tested a variant in K562 and found no allelic skew (p = 0.85), "
        "this is an empirical negative result. Allowing AlphaGenome to 'fill' this gap would overwrite ground truth with model uncertainty, waste API quota, "
        "and inflate false positive rates in downstream Bayesian fusion.",
        title="DATA INTEGRITY RULE"
    )

    # -------------------------------------------------------------
    # 6. Case Study: HBB Locus
    # -------------------------------------------------------------
    h6 = doc.add_heading("6. Case Study: The Beta-Globin (HBB) Locus Walkthrough", level=1)
    h6.paragraph_format.space_before = Pt(16)
    h6.paragraph_format.space_after = Pt(6)

    doc.add_paragraph(
        "To illustrate the behavior of Layers 0 through 2, we trace four representative variants across the human beta-globin locus (chr11:5,127,000-5,329,000, GRCh38):"
    )

    cs_tbl = doc.add_table(rows=5, cols=5)
    cs_tbl.rows[0].cells[0].paragraphs[0].text = "Variant Key"
    cs_tbl.rows[0].cells[1].paragraphs[0].text = "Biological Context"
    cs_tbl.rows[0].cells[2].paragraphs[0].text = "Layer 0 Call"
    cs_tbl.rows[0].cells[3].paragraphs[0].text = "Layer 1 Evidence"
    cs_tbl.rows[0].cells[4].paragraphs[0].text = "Layer 2 Routing & Tier"

    cs_data = [
        ("11:5227002:T:C", "HBB Promoter / 5' UTR common benign", "is_regulatory=True (promoter)", "eqtl: observed_null (MAF 0.12)\ngnomad_af: AF 0.12\nmpra/crispr: absent", "eQTL is Experimental Null -> no AG trigger; routed to Agreement Queue. MPRA/CRISPR trigger AG. Tier 2 (N_exp=1)."),
        ("11:5227172:G:A", "HBB Promoter pathogenic (beta-thalassemia)", "is_regulatory=True (5' UTR + promoter)", "eqtl: absent (MAF rare)\ngnomad_af: absent\nmpra/crispr: absent", "eQTL is Biological Gap -> triggers AG expression track. MPRA/CRISPR trigger AG. Tier 3 (N_exp=0). High-value completion candidate."),
        ("11:5280000:A:G", "Distal LCR Enhancer (50 kb upstream)", "is_regulatory=True (interval:dELS)", "mpra: observed_positive (log2FC 1.8)\neqtl: absent\ncrispr: absent", "MPRA is Positive Hit -> routed to Agreement Queue. eQTL/CRISPR trigger AG. LCR element preserved. Tier 2 (N_exp=1)."),
        ("11:5229000:G:T", "Deep Intronic / Intergenic VUS", "is_regulatory=True (interval:pELS)", "eqtl: absent\nmpra: absent\ncrispr: absent\ngwas/gnomad: absent", "Zero experimental coverage. 3 AlphaGenome triggers compiled (Whole Blood expression, K562 reporter, K562 CRISPR). Tier 3 (N_exp=0)."),
    ]

    for idx, (vk, ctx, l0, l1, l2) in enumerate(cs_data, start=1):
        cs_tbl.rows[idx].cells[0].paragraphs[0].text = vk
        cs_tbl.rows[idx].cells[1].paragraphs[0].text = ctx
        cs_tbl.rows[idx].cells[2].paragraphs[0].text = l0
        cs_tbl.rows[idx].cells[3].paragraphs[0].text = l1
        cs_tbl.rows[idx].cells[4].paragraphs[0].text = l2

    style_table(cs_tbl, col_widths=[Inches(1.2), Inches(1.3), Inches(1.3), Inches(1.3), Inches(1.4)], header_bg="1E3A8A")

    # -------------------------------------------------------------
    # 7. Engineering, Verification & Testing
    # -------------------------------------------------------------
    h7 = doc.add_heading("7. Engineering Rigor & Verification (52 Automated Tests)", level=1)
    h7.paragraph_format.space_before = Pt(16)
    h7.paragraph_format.space_after = Pt(6)

    doc.add_paragraph(
        "The codebase incorporates comprehensive automated testing using pytest. All Layer 1 and Layer 2 tests run without "
        "network access using synthetic fixtures:"
    )

    doc.add_paragraph(
        "• Layer 0 Regression Suite: Verifies LCR survival (test_lcr_survival.py), gene-restricted VEP consequences against cluster contamination (test_gene_scoped_vep.py), and alternative transcript union (test_canonical_union.py).\n"
        "• Layer 1 Unit Suite (22 tests): Verifies HTTP cache TTL and GraphQL error rejection, dynamic gnomAD rate-limiting, GWAS windowed failure isolation, GTEx pagination, and local MPRA/CRISPR schema validation.\n"
        "• Layer 2 Rule Suite (26 tests): Verifies reason categorization, multi-context status aggregation, strict exclusion of experimental nulls from AlphaGenome triggers, agreement queue segregation, and tier boundaries."
    )

    # -------------------------------------------------------------
    # 8. Roadmap & Mentor Discussion Questions
    # -------------------------------------------------------------
    h8 = doc.add_heading("8. Roadmap (Layers 3–7) & Key Questions for Mentor", level=1)
    h8.paragraph_format.space_before = Pt(16)
    h8.paragraph_format.space_after = Pt(6)

    doc.add_paragraph(
        "Next Development Phases:\n"
        "• Layer 3 (In Progress): AlphaGenome API client binding, 1 Mb receptive field sequence extraction, batch inference execution.\n"
        "• Layer 4 & 5: Delta-track normalization and Hierarchical Bayesian Evidence Fusion.\n"
        "• Layer 6: ACMG-style calibration across Tiers 1-3 using ClinVar benchmarks.\n"
        "• Layer 7: Clinical prioritization reporting and audit manifests."
    )

    doc.add_heading("Strategic Discussion Questions for Our Mentor:", level=2)
    doc.add_paragraph(
        "1. AlphaGenome Proxy Penalization: When using AlphaGenome chromatin accessibility as a surrogate for episomal MPRA, "
        "should we use a fixed Dirichlet penalty factor or learn empirical shrinkage weights from the Layer 2 Agreement Queue?\n"
        "2. Prior Specification in Evidence-Sparse Tier 3: In Tier 3 (zero experimental coverage), the baseline clinical prior is low (~0.001). "
        "Should this prior be uniform or modulated by ENCODE cCRE functional categories (promoter vs distal enhancer)?\n"
        "3. Pleiotropic Tissue Aggregation: For genes expressed in multiple GTEx tissues, should Bayesian fusion evaluate the maximum tissue effect "
        "or a tissue-covariance-weighted composite?\n"
        "4. Secondary Benchmark Loci: Beyond HBB, which disease loci (e.g. LDLR, TERT, CFTR, BRCA1) would you recommend prioritizing for cross-locus validation?"
    )

    # Save document
    out_dir = "docs"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "mentor_presentation_layers_0_to_2.docx")
    doc.save(out_path)
    print(f"Successfully generated presentation Word document: {out_path}")
    return out_path

if __name__ == "__main__":
    build_word_document()
