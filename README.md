# GBM WEBS Project - Reproducible Analysis Pipeline

**Paper**: Post-transcriptional spliceosome accumulation marks a writer-dominant malignant cell state in glioblastoma

**Authors**: Jia Chen, Junhan Yang, Bing Yang, Tianyi Zhou, Peijun Zhao

**Corresponding author**: Jia Chen, Department of Bioinformatics, Chengdu Neusoft University, Chengdu, China (1176625889@qq.com)

---

## Overview

This repository contains the complete analysis pipeline for the manuscript. All analyses were re-executed from raw public data on an independent machine (R 4.5.2).

**Core finding**: A writer-dominant expression state (quantified by the Writer–Eraser Balance Score, WEBS) is associated with proteome-wide spliceosome accumulation in GBM (107/111 proteins up-regulated; protein-layer ρ = 0.79 vs mRNA-layer ρ = −0.03), accompanied by coordinated acetylation/phosphorylation and altered alternative splicing output (4,387 events, replicated in TCGA-LGG).

## Project Structure

```
GBM_WEBS_Project/
├── 00_setup.R                    # Install all R packages
├── 01_download_data.R            # Download public datasets
├── 02_scrna_analysis.R           # Single-cell analysis (core pipeline)
├── 03_pseudotime.R               # Slingshot trajectory analysis
├── 04_spatial.R                  # 10x Visium spatial transcriptomics
├── 05_snatac.R                   # snATAC-seq analysis (Signac)
├── 06_pathway_tf.R               # AUCell pathway + SCENIC TF analysis
├── 07_bulk_survival.R            # Cohort survival analysis
├── 08_cptac_proteomics.R         # CPTAC TMT proteomics
├── 09_sixgene_signature.R        # LASSO-Cox risk signature
├── 10_pharmacogenomics.R         # DepMap + GDSC2 screening
├── 11_benchmarking.R             # Method comparison (ssGSEA/AUCell/PLAGE)
├── 12_figures.R                  # Standardized plotting functions
├── 13-65_*.R                     # Extended analyses (see below)
├── tools/
│   └── WEBS/
│       └── webs.R                # Standalone WEBS implementation
└── results/ figures_final/       # Outputs
```

## Extended analysis scripts

| Script | Analysis |
|--------|----------|
| 14 | Parameter sensitivity (maxRank, leave-one/two-out) |
| 15 | Bootstrap stability of feature selection (1,000 iterations) |
| 16 | Computational performance benchmark |
| 17 | STRING PPI network |
| 18 | Spatial autocorrelation (Moran's I) |
| 19 | SHAP interpretability |
| 20 | Acetylation site-level analysis |
| 21 | CopyKAT purity control |
| 22 | CellChat communication |
| 23 | MOFA2 multi-omics integration |
| 27 | RPPA lactate-metabolism score |
| 28b | Specificity control (10,000 expression-matched random gene sets) |
| 38 | Confounder adjustment, subgroups, treatment interaction |
| 40/51 | Pan-cancer analysis (33 TCGA cohorts) |
| 44 | Virtual knockout (GENIE3 GRN) |
| 45 | GSE83300 external validation |
| 47 | GLASS validation + longitudinal pairs |
| 48 | Alternative splicing (TCGASpliceSeq PSI vs WEBS) |
| 49 | Phosphoproteome analysis |
| 50 | GLASS recurrence-interval association |
| 55 | Pseudoreplication control (pseudobulk + mixed-effects) |
| 56 | Splicing immune-fraction correction |
| 57 | PTM parent-protein-abundance correction |
| 58 | RNA–protein uncoupling rigor (bootstrap CI, LOO) |
| 59 | WEBS incremental value beyond cell-state programs |
| 60 | Cross-cohort splicing replication (TCGA-LGG) |
| 61 | Signature performance (C-index, time-dependent AUC) |
| 62 | 107/111 threshold sensitivity |
| 65 | GLASS first-pair-per-patient sensitivity |

## Data sources

All datasets are public: GSE138794, GSE194329, GSE131928, GSE7696, GSE83300 (GEO); TCGA cohorts (GDC / UCSC Xena GDC hub); CGGA-693/325 (cgga.org.cn, registration required); CPTAC-GBM (LinkedOmics); TCGASpliceSeq; DepMap 24Q2; GDSC2. GLASS data were accessed via Synapse (syn17038081) under the consortium's data-use terms and are NOT redistributed here.

## WEBS standalone tool

`tools/WEBS/webs.R` provides the self-contained WEBS implementation for single-cell and bulk expression data.

## Citation

If you use WEBS or this pipeline, please cite the manuscript (under review) and archive DOI: Zenodo (upon publication)
