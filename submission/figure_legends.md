# Figure Legends

**Figure 1 | A computational framework for discovering a therapeutically targetable malignant cell state in glioblastoma.** Six data sources (scRNA-seq of 35,823 cells, spatial transcriptomics, snATAC-seq, eight bulk cohorts, CPTAC proteomics, DepMap/GDSC2) feed the Writer–Eraser Balance Score (WEBS), followed by a three-tier analysis (Discovery: cell-state landscape and hypoxia decoupling; Mechanism: post-transcriptional spliceosome accumulation and acetylation coupling; Translation: six-gene signature, pharmacogenomics, molecular docking). Key finding callout: protein ρ = +0.79 vs mRNA ρ = −0.03.

**Figure 2 | The Writer–Eraser Balance Score (WEBS) landscape.** (a) Cell-type annotation of 35,823 cells from 20 GBM tumor samples (marker-based with original-study labels). (b) WEBS projected onto the UMAP embedding. (c) WEBS distribution across the IQR-defined groups (LLS, bottom 25%; LDTS, middle 50%; LHS, top 25%; Kruskal–Wallis p < 10⁻³⁰⁰).

**Figure 3 | Neftel states and WEBS.** (A) Neftel module scores (MES/AC/NPC/OPC) by WEBS group in 7,952 neoplastic cells. (B) WEBS vs mesenchymal module score with GAM fit. External validation in GSE131928 (Smart-seq2: Kruskal–Wallis p = 1.6 × 10⁻¹⁶⁸; 10x: p = 1.0 × 10⁻³¹³).

**Figure 4 | WEBS is decoupled from canonical hypoxia markers.** (a) WEBS across invasion region types (all cells; Kruskal–Wallis p < 2.2 × 10⁻¹⁶). (b) Patient-level paired Margin vs Core WEBS (paired Wilcoxon, n = 17). (c) WEBS vs CA9 + VEGFA score (all cells, n = 35,823; within malignant cells ρ = −0.054). (d) Hallmark Hypoxia score by WEBS group; the Margin-like association persists after hypoxia adjustment (logistic regression p = 3.8 × 10⁻⁴).

**Figure 5 | Pseudotime trajectory connects LLS to LHS.** Slingshot trajectory on neoplastic cells (LLS origin): UMAP colored by pseudotime (left), WEBS (center), and GAM fit of WEBS along pseudotime (right; ρ = 0.301, p = 2.7 × 10⁻¹⁶⁶).

**Figure 6 | Cell–cell communication across WEBS states.** (A) Differential interaction number and strength between LLS and LHS (CellChat v2). (B) Pathway-level information-flow comparison; SPP1–CD44 up in LHS, MIF–CD74/CXCR4 enriched in LLS (permutation p < 10⁻⁴). (C) Differential interaction network (red: stronger in LHS).

**Figure 7 | Spatial transcriptomics.** Per-spot WEBS on six 10x Visium GBM slices (symmetric blue–white–red scale centered at 0). WEBS hotspots show significant spatial autocorrelation in 3 of 6 slices (Moran's I, permutation p ≈ 0.001).

**Figure 8 | snATAC-seq promoter accessibility.** Accessibility at writer/eraser promoter peaks (TSS ± 2 kb) across six snATAC libraries from GSE138794. (Fragment files were not deposited for this dataset; analysis is limited to peak-matrix level.)

**Figure 9 | Functional characterization of WEBS states.** (A) Hallmark pathway activity differences (AUCell, LHS − LLS) in malignant cells. (B) Differential transcription-factor expression between LHS and LLS.

**Figure 10 | External validation and six-gene signature.** (A) WEBS by Neftel state in GSE131928. (B) Cross-platform patient-level consistency (Smart-seq2 vs 10x). (C–E) Kaplan–Meier curves for the six-gene signature in TCGA-GBM (p = 0.0028), CGGA-693 (p = 3.4 × 10⁻¹⁹), and CGGA-325 (p = 7.7 × 10⁻¹⁸).

**Figure 11 | Pharmacogenomic prioritization.** (A) WEBS vs CRISPR dependency for MED1, SF3B1, MLLT1, KAT2A (DepMap 24Q2, 102 GBM/glioma lines); MLLT1/ENL did not show GBM-selective essentiality and was not prioritized for therapeutic follow-up. (B) HDAC-inhibitor sensitivity vs WEBS (GDSC2 vorinostat/entinostat). (C) Top WEBS-sensitizing compounds. (D) Molecular docking: vorinostat–HDAC2 (−7.9 kcal/mol), H3B-8800–SF3B1 (−6.7 kcal/mol), AGK2–SIRT2 (−8.5 kcal/mol).

**Figure 12 | Working model.** Writer–eraser imbalance drives splicing-factor acetylation, post-transcriptional spliceosome accumulation, and therapeutic vulnerabilities. (Original figure retained from v16; update to include the acetylation mechanism is recommended.)

---

**Table 1 | Cohorts used in this study.** Platform, sample size, and analytical role for all 15 datasets.

**Table 2 | Therapeutic target prioritization.** Tiered ranking by convergent evidence (proteomics, CRISPR dependency, drug sensitivity, bootstrap stability, docking).
