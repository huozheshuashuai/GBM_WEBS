# GBM WEBS Project - Reproducible Analysis Pipeline

**Paper**: Spliceosome activation defines an epigenetically hyperactive and therapeutically targetable malignant cell state in glioblastoma

**Authors**: Jia Chen, Junhan Yang, Bing Yang, Tianyi Zhou, Lu Zhang

**Institution**: Chengdu Neusoft University

---

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
├── 07_bulk_survival.R            # 8-cohort survival analysis
├── 08_cptac_proteomics.R         # CPTAC TMT proteomics
├── 09_sixgene_signature.R        # LASSO-Cox risk signature
├── 10_pharmacogenomics.R         # DepMap + GDSC2 screening
├── 11_benchmarking.R             # Method comparison (ssGSEA/AUCell/PLAGE)
├── 12_figures.R                  # Standardized plotting functions
├── tools/
│   └── WEBS/
│       └── webs.R                # Standalone WEBS implementation
├── data/                         # Data directory (created by 01_download_data.R)
│   ├── raw/
│   └── processed/
├── results/                      # Analysis outputs
└── figures/                      # Generated figures
```

---

## Quick Start

### Step 1: Install Dependencies

```r
source("00_setup.R")
```

Requires **R >= 4.3.0**, preferably **R 4.5.0**.

### Step 2: Download Data

```r
source("01_download_data.R")
# Run individual download functions as needed:
# download_gse138794()
# download_gse194329()
# download_tcga_gbm()
# etc.
```

**Note**: Some datasets require manual download due to size limitations:
- GSE138794 count matrices (~GB scale)
- CPTAC proteomics (TMT data)
- DepMap CRISPR files

### Step 3: Run Analysis Pipeline

```r
# Load scRNA-seq data (replace with your actual data path)
obj <- load_scrna("data/raw/GSE138794")

# Run full pipeline: QC -> Normalize -> Cluster -> Annotate -> WEBS
obj <- run_full_pipeline(obj)

# Save
saveRDS(obj, "data/processed/GSE138794_webs_annotated.rds")
```

### Step 4: Run Individual Modules

```r
# Pseudotime analysis
source("03_pseudotime.R")
pt_results <- run_slingshot(obj)

# Spatial analysis
source("04_spatial.R")
spatial_obj <- load_visium("data/raw/GSE194329/sample1", "sample1")

# Survival analysis
source("07_bulk_survival.R")
clinical <- compute_bulk_webs(expr_matrix, clinical_data)
plot_km(clinical, "TCGA-GBM")
```

---

## WEBS Usage (Standalone)

```r
source("tools/WEBS/webs.R")

# Compute WEBS on Seurat object
obj <- webs_compute(obj, 
                    writer_genes = c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP"),
                    eraser_genes = c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2",
                                     "SIRT3", "SIRT6", "SIRT7"))

# Classify cells
obj <- webs_classify(obj, method = "iqr")

# Plot distribution
webs_plot(obj, group.by = "cell_type")
```

---

## Datasets Summary

| Dataset | Source | N | Platform | Usage |
|---------|--------|---|----------|-------|
| GSE138794 | GEO | 35,635 cells | 10x Genomics | Discovery scRNA-seq |
| GSE194329 | GEO | 6 slices | 10x Visium | Spatial transcriptomics |
| GSE131928 | GEO | 24,131 cells | Smart-seq2 + 10x | External validation |
| TCGA-GBM | cBioPortal | ~150 samples | RNA-seq | Survival |
| CGGA-693 | CGGA | 693 samples | RNA-seq | LASSO training |
| CGGA-156 | CGGA | 156 samples | RNA-seq | Validation |
| GSE7696 | GEO | ~80 samples | Microarray | Survival |
| GSE13041 | GEO | ~200 samples | Microarray | Survival |
| REMBRANDT | GEO | ~200 samples | Microarray | Survival |
| GLASS | GLASS consortium | ~200 samples | RNA-seq | Recurrence |
| CPTAC | LinkedOmics | 99 samples | TMT proteomics | Protein-level |
| DepMap 24Q2 | DepMap | 1000+ lines | CRISPR | Dependency |
| GDSC2 | Sanger | 1000+ lines | Drug screen | Sensitivity |
| IMvigor210 | easierData | ~300 samples | RNA-seq | Immunotherapy |

---

## Key Parameters

| Analysis | Key Parameters | R Package |
|----------|---------------|-----------|
| WEBS | UCell rank-based AUC | UCell v2.6.2 |
| scRNA QC | >200 genes, >500 UMIs, <20% MT | Seurat v5.1.0 |
| Clustering | PCA (1:30), resolution=0.8 | Seurat v5.1.0 |
| Pseudotime | UMAP-based, Core=start | Slingshot v2.14.0 |
| Pathways | leading-edge=0.375 | AUCell v1.26.0 |
| Differential | Wilcoxon (sc) / Limma (bulk) | Seurat / limma |
| Survival | Cox PH, Harrell's C-index | survival v3.7-0 |
| Signature | LASSO-Cox, alpha=1, 10-fold CV | glmnet v4.1-8 |
| Protein DE | Limma-trend | limma v3.62.0 |

---

## Software Versions

- R: 4.5.0
- Seurat: 5.1.0
- SeuratObject: 5.0.2
- Signac: 1.14.0
- UCell: 2.6.2
- glmnet: 4.1-8
- survival: 3.7-0
- limma: 3.62.0
- AUCell: 1.26.0

---

## Citation

If you use this code, please cite:

> Chen J, Yang J, Yang B, Zhou T, Zhang L. Spliceosome activation defines an epigenetically hyperactive and therapeutically targetable malignant cell state in glioblastoma. *BBA - Molecular Basis of Disease* (submitted).

---

## Contact

Corresponding author: Jia Chen (1176625889@qq.com)

Chengdu Neusoft University, Department of Health Science and Medical Engineering
