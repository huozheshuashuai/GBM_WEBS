# =============================================================================
# Data Download Scripts
# =============================================================================
# All datasets are publicly available. Run each section as needed.
# =============================================================================

library(GEOquery)
library(TCGAbiolinks)

# ---- 1. scRNA-seq & snATAC-seq (Discovery Cohort) ----
# GSE138794: 35,635 cells from 17 GBM patients (19 tumor samples)
# Paper: Neftel et al. (Cell, 2019) + additional samples
# Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138794

download_gse138794 <- function(outdir = "data/raw/GSE138794") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Download GEO series matrix and supplementary files
  gse <- getGEO("GSE138794", GSEMatrix = TRUE, getGPL = FALSE)
  
  # The count matrices are typically in supplementary files
  # You may need to manually download from GEO
  message("GSE138794 metadata downloaded. Count matrices require manual download from:")
  message("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138794")
  message("Look for 'GSE138794_GBM_all_cells_10x_counts.RDS' or similar")
  
  return(gse)
}

# ---- 2. Spatial Transcriptomics ----
# GSE194329: Six 10x Genomics Visium slices
# Paper: Ren et al. (Cancer Cell, 2023)
# Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194329

download_gse194329 <- function(outdir = "data/raw/GSE194329") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  gse <- getGEO("GSE194329", GSEMatrix = TRUE, getGPL = FALSE)
  message("GSE194329 metadata downloaded. Spatial data (H5 files) require manual download.")
  return(gse)
}

# ---- 3. External scRNA-seq Validation ----
# GSE131928: Neftel et al. (Cell, 2019) - 24,131 cells from 28 patients
# Smart-seq2 + 10x Genomics
# Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928

download_gse131928 <- function(outdir = "data/raw/GSE131928") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  gse <- getGEO("GSE131928", GSEMatrix = TRUE, getGPL = FALSE)
  message("GSE131928 metadata downloaded. Expression matrices require manual download.")
  return(gse)
}

# ---- 4. TCGA-GBM Bulk RNA-seq ----
# Via cBioPortal or TCGAbiolinks
# cBioPortal: https://www.cbioportal.org/study/summary?id=gbm_tcga

download_tcga_gbm <- function(outdir = "data/raw/TCGA_GBM") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  query <- GDCquery(
    project = "TCGA-GBM",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  GDCdownload(query, directory = outdir)
  
  data <- GDCprepare(query, directory = outdir)
  saveRDS(data, file.path(outdir, "TCGA_GBM_RNAseq.rds"))
  return(data)
}

# ---- 5. CGGA Datasets ----
# CGGA-693 and CGGA-156
# Download from: http://www.cgga.org.cn/
# Requires registration

# ---- 6. Other GEO Cohorts ----
# GSE7696, GSE13041, REMBRANDT
# Download individually via GEOquery

download_geo_cohort <- function(accession, outdir = NULL) {
  if (is.null(outdir)) outdir <- paste0("data/raw/", accession)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  gse <- getGEO(accession, GSEMatrix = TRUE, getGPL = FALSE)
  saveRDS(gse, file.path(outdir, paste0(accession, "_metadata.rds")))
  message(accession, " metadata downloaded.")
  return(gse)
}

# ---- 7. CPTAC Proteomics ----
# LinkedOmics: http://linkedomics.org/data_download/TCGA-GBM/
# Or via Broad Institute CPTAC portal
# Files needed:
#   - Gene-level proteome (11,141 proteins x 99 samples)
#   - Site-level acetylome (13,556 sites x 99 samples)

download_cptac_gbm <- function(outdir = "data/raw/CPTAC") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("CPTAC data requires manual download from:")
  message("https://proteomics.cancer.gov/data-portal")
  message("Or: http://linkedomics.org/data_download/TCGA-GBM/")
  message("Download GBM proteome and acetylome TMT datasets")
}

# ---- 8. DepMap 24Q2 ----
# https://depmap.org/portal/download/all/
# Files needed:
#   - CRISPR_gene_effect.csv (Chronos scores)
#   - Model.csv (cell line metadata)
#   - OmicsExpressionProteinCodingGenesTPMLogp1.csv (expression)

download_depmap <- function(outdir = "data/raw/DepMap") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("DepMap data requires manual download from:")
  message("https://depmap.org/portal/download/all/")
  message("Download 24Q2 release: CRISPRGeneEffect, Model, OmicsExpression")
}

# ---- 9. GDSC2 ----
# https://www.cancerrxgene.org/downloads/bulk_download
# Files needed:
#   - GDSC2_fitted_dose_response_17Jul24.txt
#   - Cell_line_RMA_proc_basalExp.txt

download_gdsc2 <- function(outdir = "data/raw/GDSC2") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("GDSC2 data requires manual download from:")
  message("https://www.cancerrxgene.org/downloads/bulk_download")
  message("Download GDSC2 fitted dose response and expression data")
}

# ---- 10. IMvigor210 (Immunotherapy) ----
# Package: easierData / ExperimentHub
# Mariathasan et al. (Nature, 2018)

download_imvigor210 <- function() {
  if (!require("easierData", quietly = TRUE)) {
    BiocManager::install("easierData")
  }
  library(ExperimentHub)
  eh <- ExperimentHub()
  query(eh, "Mariathasan2018_PDL1_treatment")
  message("IMvigor210 data accessed via easierData/ExperimentHub")
}

# ---- RUN ALL ----
# Uncomment lines below to download
# download_gse138794()
# download_gse194329()
# download_gse131928()
# download_tcga_gbm()
# download_geo_cohort("GSE7696")
# download_geo_cohort("GSE13041")
# download_cptac_gbm()
# download_depmap()
# download_gdsc2()
# download_imvigor210()
