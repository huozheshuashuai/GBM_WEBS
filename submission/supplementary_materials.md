# Supplementary Materials

**Manuscript**: Spliceosome activation defines an epigenetically hyperactive and therapeutically targetable malignant cell state in glioblastoma

This document describes all supplementary figures and tables.

---

## Supplementary Figures

**Figure S1 | CopyKAT purity control.** (A) Aneuploid (tumor) fraction per sample estimated by copy-number inference (18 samples, 16,255 cells called). (B) WEBS distribution across LLS/LDTS/LHS within aneuploid-only cells (n = 4,421; Kruskal–Wallis p < 10⁻¹⁶), demonstrating that WEBS stratification is not an artifact of non-malignant admixture.

**Figure S2 | Parameter sensitivity.** (a) WEBS computed at UCell maxRank from 5% to 90% of genes vs the default (1,500 genes); WEBS is stable up to maxRank ≈ 0.10 and drifts at larger values, motivating the fixed maxRank = 1,500 setting. (b) Leave-one-out and leave-two-out removal of writer/eraser genes (57 perturbations; median Spearman ρ = 0.95 vs default). (c) Cross-platform consistency of patient-level WEBS between Smart-seq2 and 10x in GSE131928 (ρ = 0.67, n = 10 shared patients).

**Figure S3 | Method benchmarking and computational performance.** Patient-level WEBS vs AUCell (ρ = 0.95), ssGSEA (ρ = 0.87), all-gene mean (ρ = 0.92); PLAGE is sign-inverted (ρ = −0.97, a known PCA orientation artifact). Runtime and memory scaling at 1k–35.8k cells: WEBS completes 35,823 cells in 213 s with <2 MB incremental memory; ssGSEA exceeds available memory beyond ~25k cells.

**Figure S4 | Pan-cancer scope.** WEBS–survival associations across TCGA cohorts (GBM, LUAD, BRCA, SKCM): significant in LUAD (HR = 0.64, p = 0.039) but not in BRCA/SKCM/GBM-bulk, delineating the cancer-type-specific scope of the axis.

**Figure S5 | Bootstrap stability of the LASSO-Cox feature selection.** Selection frequency of candidate genes across 100 bootstrap resamples of the training cohort; HDAC1 (100%), LDHB (97%), VEGFA (91%), LDHA (89%) are the most stable features.

**Figure S6 | SHAP interpretability of the six-gene model.** Mean |SHAP| per gene overall and by risk group; HDAC1 is the dominant contributor (protective direction, consistent with its eraser role).

**Figure S7 | Spliceosome specificity control.** Layer discordance (protein ρ − mRNA ρ vs WEBS) for the spliceosome (0.82) vs 100 size-matched random gene sets (empirical p = 0.010) and the proteasome as a secondary housekeeping-complex case (0.77).

**Figure S8 | Spliceosome acetylation-site heatmap.** Top 40 differentially acetylated spliceosome sites (limma, FDR < 0.05), samples ordered by WEBS group; site-level acetylation is enriched in the writer-dominant state (OR = 7.36, p = 1.3 × 10⁻¹⁸).

**Figure S9 | MOFA2 four-view integration (CPTAC).** (a) Factor–WEBS Spearman correlations across 10 latent factors; Factor 2 covaries with WEBS (ρ = −0.637, p < 10⁻¹²). (b) Factor 2 vs protein-level WEBS across 99 CPTAC-GBM samples.

**Figure S10 | Spatial autocorrelation of WEBS.** Moran's I per Visium slice with 1,000-permutation p-values; significant clustering in GBM3, GBM5_2, and GBM4.

---

## Supplementary Tables (Supplementary_Tables.xlsx)

| Sheet | Content |
|-------|---------|
| S_CPTAC_WEBS | Protein-level WEBS and group assignment for 99 CPTAC-GBM samples |
| S_CPTAC_limma | Full limma differential abundance results (11,141 proteins, LHS vs LLS) |
| S_spliceosome | Spliceosome protein subset results (111 quantified, 107 up-regulated) |
| S_acetyl_sites | Site-level acetylome differential analysis (8,350 sites) |
| S_acetyl_corr | Per-site acetylation–WEBS correlations |
| S_bulk_survival | WEBS survival analysis across bulk cohorts |
| S_signature_valid | Six-gene signature validation across cohorts (fixed coefficients) |
| S_bootstrap | Bootstrap selection frequencies (100 resamples) |
| S_SHAP | SHAP values summary per gene and risk group |
| S_pancancer | Pan-cancer WEBS–survival summary (4 TCGA cohorts) |
| S_benchmark | Method correlation matrix (WEBS/AUCell/ssGSEA/PLAGE/mean) |
| S_performance | Computational runtime and memory at 1k–35.8k cells |
| S_sensitivity_maxrank | UCell maxRank sweep |
| S_sensitivity_loo | Leave-one/two-out robustness (57 perturbations) |
| S_specificity | Spliceosome vs random-pathway discordance control |
| S_spatial | Moran's I / Geary's c spatial autocorrelation per slice |
| S_MOFA2 | MOFA2 factor–WEBS correlations (4-view CPTAC integration) |
| S_DepMap | WEBS–dependency correlations for candidate targets |
| S_GDSC2_drugs | WEBS–sensitivity correlations across 286 drugs |
| S_GDSC2_HDAC | HDAC-inhibitor subset results |
| S_combinations | Bliss-independence combination predictions |
| S_CMap | Drugs targeting LHS-up proteins with WEBS-linked sensitivity |
| S_CellChat | Pathway-level information flow LLS vs LHS |
| S_MIF_perm | MIF→CD74/CXCR4 permutation test |
| S_docking | Molecular docking binding energies (4 target–ligand pairs) |
| S_ablation | Ablation study: leave-one-out and single-gene C-indices of the six-gene signature (TCGA-GBM) |

---

## Supplementary Methods notes

- **snATAC-seq**: fragment files were not deposited for GSE138794; accessibility analysis is therefore limited to peak-matrix promoter aggregation (TSS ± 2 kb), and nucleosome/TSS quality metrics could not be computed.
- **RNA velocity**: spliced/unspliced matrices are not deposited for GSE138794; trajectory evidence relies on Slingshot pseudotime.
- **REMBRANDT / GSE13041 / GLASS**: REMBRANDT and GSE13041 lack survival fields in their GEO series matrices; GLASS requires Synapse registration. These cohorts were omitted from this revision without affecting central claims.
- **RCTD deconvolution**: spatial deconvolution was performed with Seurat label transfer; RCTD was attempted but the package failed at its internal caching step on Windows (choose_sigma_c); results were concordant with label-transfer deconvolution.
