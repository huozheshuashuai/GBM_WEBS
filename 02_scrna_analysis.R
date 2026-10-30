# =============================================================================
# Single-Cell RNA-seq Analysis Pipeline (Discovery Cohort: GSE138794)
# =============================================================================
# Paper: 35,635 cells from 17 GBM patients (19 tumor samples)
# Platform: 10x Genomics Chromium
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ---- Configuration ----
WRITER_GENES <- c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP")
ERASER_GENES <- c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")

CELL_TYPE_MARKERS <- list(
  Neuron = c("RBFOX3", "MAP2"),
  Astrocyte = c("GFAP"),
  Macrophage_Microglia = c("PTPRC", "CD14", "CSF1R", "C1QA", "AIF1"),
  Endothelial = c("PECAM1"),
  Pericyte = c("PDGFRB"),
  Oligodendrocyte = c("OLIG1", "OLIG2", "MBP"),
  Malignant = c("SOX2", "EGFR", "OLIG2")
)

NEFTEL_MODULES <- list(
  MES1 = c("CD44", "VIM", "S100A11", "S100A10", "ANXA1", "LGALS1", "MT2A", "C1S", "NAMPT", "EMP3"),
  MES2 = c("BCL2A1", "CD83", "CXCL8", "AREG", "CD74", "TIMP1", "CCL2", "IL1B", "CCL4", "C3"),
  AC = c("APOE", "GFAP", "S100B", "HEPN1", "ARHGEF26", "FAM107A", "GJA1", "SPARCL1", "BMPR1B", "HP", "GLIPR1", "CLU", "AGT", "PMP2", "ABCA1", "DBI", "EDNRB", "CLDN10", "OLR1", "PLTP", "FABP7", "GATM", "HOPX", "MLC1", "PLP1", "CST3", "TIMP4", "ALDOC", "NTRK2", "RAMP1", "PGCP", "PTPRZ1", "TTYH1", "SERPINF1", "ID4", "NDRG2", "CST1", "SLC1A3", "ITPR2", "EMP1", "ATP1A2", "CRYAB", "RRBP1", "SLC14A1", "ADM", "PANTR1", "CD81", "EGFR", "BNIP3", "AHNAK2", "GRN", "PLOD1"),
  NPC1 = c("SOX4", "STMN1", "STMN2", "BEX1", "SOX11", "HMP19", "TMSB10", "DLL3", "MARCKSL1", "PTPRS", "GAP43", "BEX4", "INA", "NNAT", "TUBB2B", "DCX", "CD24", "STMN4", "SNAP25", "NPPA", "RBFOX2", "tagln2", "MAP1B", "MDK", "TCIRG1", "PGAM1", "PCP4", "CDK4", "NMB", "MIAT", "NRCAM", "LGALS3", "BASP1", "CRABP2", "SNCG", "YWHAE", "TTYH2", "TPT1", "RTN1", "GNB1", "PPP1R14B", "ENO2", "TUBA1A", "NSG1", "DBI", "PFN2", "INA", "PKM", "TUBB3", "NAP1L5", "NRN1", "HNRNPA3", "ELAVL4", "ATP6V1C1", "TCEAL3", "NANOGP8", "RAB3A"),
  NPC2 = c("VGF", "SOX11", "DCX", "DLL3", "PTPRS", "SNAP25", "STMN1", "BASP1", "CD24", "INA", "STMN2", "MAP1B", "SOX4", "GAP43", "GNG3", "HMP19", "PFN2", "STMN4", "TAGLN3", "TUBB2B", "DBI", "CRABP2", "BEX1", "MARCKSL1", "NRCAM", "LGALS3", "NAP1L5", "GNB1", "TMSB10", "ATPIF1", "MDK", "RTN1", "ENO2", "PPP1R14B", "TCEAL3", "SNCG", "TCIRG1", "ATP6V1C1", "RAB3A", "NANOGP8", "NSG1", "TPT1", "BEX4", "NMB", "PGAM1", "NRN1", "TUBB3", "CDK4", "ELAVL4", "TTYH2", "YWHAE", "PCP4", "ALKBH5", "HNRNPA3", "NREP", "PKM"),
  OPC = c("SOX10", "PLP1", "BCAS1", "S100B", "GPR17", "TF", "PTPRZ1", "CA2", "DBI", "MYRF", "MAG", "MOG", "APOD", "ERMN", "CLDN11", "ANLN", "KLK6", "TMEM144", "MBP", "PLLP", "FA2H", "MOBP", "CRYAB", "TRF", "SIRT2", "RALA", "RTN4", "TPPP3", "TSPAN2", "ABCA2", "OPALIN", "CAR2", "QDPR", "AATK", "CNP", "GJB1", "TFEB", "ST18", "PLEKHB1", "CERCAM", "SLC44A1", "RASSF2", "BIN1", "PEX5L", "ALCAM", "SPARCL1", "CASPR", "SPOCK3", "NKAIN4", "NKX6-2", "TEX14", "PALM2", "LGI3", "PPP1R16B", "TMEM88B", "PADI2", "SEPP1", "SERINC5", "LAMP1", "TMEFF2", "CLDN11", "BMPR1B", "RAB31", "SNX1", "CYB5R2", "HEPACAM", "STXBP3", "S100A16", "LPAR1", "MAN2A1", "TPR", "EDIL3", "MAL", "CWH43", "SLC12A2", "RHOA", "EPS8", "TCN2", "CNP", "S100A13", "MARCKSL1", "CSPG5", "GPR37L1", "CRYM", "LIPA", "TBCB", "CMTM5", "TMEM125", "RFFL", "KCNJ10", "SH3GL3", "FA2H", "TPPP3", "APLP1", "EDNRB", "SPARCL1", "PPP1R14B", "ALDH1A1", "MARCKSL1", "SLC1A3", "ADAM22", "GRIA2", "CLU", "EMP2", "SELENOP", "PTGDS", "ATP1B2", "PTPRZ1", "SCD5", "PRR18", "FTH1", "BMPR1B", "CA2", "VOPP1", "APOD", "BCAS1", "TIMP4", "SOX6", "CRYAB", "HEPN1", "PLP1", "AGT", "PDGFRA", "OLIG1", "SOX8", "CSPG4", "TCN2", "ABCG2", "GPR17", "THBS3", "FGFR3", "FAM107A", "RGR", "C1QL1", "GJB1", "GPM6B", "SLC44A1", "ARL4C", "SLC1A1", "EPB41L2", "ERMN", "ATP1A2", "PTPRZ1", "TTYH2", "BCAN", "MAG", "TRF", "CA2", "APOD", "LPAR1", "MYRF", "MAL", "PLLP", "MOG", "S100B", "DBI", "FA2H", "OPALIN", "MOBP", "QDPR", "CRYAB", "ANK3", "CNP", "TMEM125", "SPARCL1", "CLDN11", "ANLN", "TMEM144", "KLK6", "SIRT2", "RALA", "RTN4", "TSPAN2", "ABCA2", "CERCAM", "BIN1", "NKAIN4", "NKX6-2", "TEX14", "LGI3", "PADI2", "SEPP1", "SERINC5", "LAMP1", "ST18", "PLEKHB1", "RASSF2", "PEX5L", "ALCAM", "CASPR", "SPOCK3", "PPP1R16B", "TMEM88B", "HEPACAM", "STXBP3", "S100A16", "MAN2A1", "TPR", "EDIL3", "CWH43", "SLC12A2", "RHOA", "EPS8", "S100A13", "CSPG5", "GPR37L1", "CRYM", "LIPA", "TBCB", "CMTM5", "KCNJ10", "SH3GL3", "APLP1", "ALDH1A1", "ADAM22", "GRIA2", "EMP2", "SELENOP", "PTGDS", "ATP1B2", "SCD5", "PRR18", "FTH1", "VOPP1", "SOX6", "HEPN1", "PDGFRA", "SOX8", "CSPG4", "THBS3", "FGFR3", "RGR", "C1QL1", "GPM6B", "ARL4C", "SLC1A1", "EPB41L2", "BCAN")
)

# ---- 1. Load Data ----
load_scrna <- function(path = "data/raw/GSE138794") {
  # Expects Seurat object or 10x data
  # If starting from 10x output: CreateSeuratObject(counts = Read10X(...))
  message("Loading GSE138794 scRNA-seq data...")
  
  # Example: if you have a pre-built RDS
  # obj <- readRDS(file.path(path, "GSE138794_seurat.rds"))
  
  # Example: from 10x raw counts
  # counts <- Read10X(data.dir = path)
  # obj <- CreateSeuratObject(counts = counts, project = "GBM", min.cells = 3, min.features = 200)
  
  message("Please load your data here based on format (RDS/10X/H5AD)")
  return(NULL)
}

# ---- 2. Quality Control ----
run_qc <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  
  # QC thresholds from Methods:
  # >200 detected genes, >500 UMIs, <20% mitochondrial reads
  obj <- subset(obj, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 20)
  
  message("After QC: ", ncol(obj), " cells retained")
  return(obj)
}

# ---- 3. Normalization & Scaling ----
preprocess <- function(obj, nfeatures = 2000) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  obj <- ScaleData(obj)
  return(obj)
}

# ---- 4. Dimensionality Reduction & Clustering ----
run_dimred <- function(obj, npcs = 30, resolution = 0.8) {
  obj <- RunPCA(obj, features = VariableFeatures(object = obj))
  obj <- FindNeighbors(obj, dims = 1:npcs)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj, dims = 1:npcs)
  return(obj)
}

# ---- 5. Cell Type Annotation ----
annotate_cells <- function(obj) {
  # Score each cell type
  for (ct in names(CELL_TYPE_MARKERS)) {
    markers <- intersect(CELL_TYPE_MARKERS[[ct]], rownames(obj))
    if (length(markers) > 0) {
      obj <- AddModuleScore(obj, features = list(markers), name = ct)
      # Rename column
      colnames(obj@meta.data)[colnames(obj@meta.data) == paste0(ct, "1")] <- ct
    }
  }
  
  # Simple annotation: assign to max score
  score_cols <- names(CELL_TYPE_MARKERS)
  scores <- obj@meta.data[, score_cols, drop = FALSE]
  obj$cell_type <- colnames(scores)[max.col(scores)]
  
  # Refine: if Malignant score is highest, call it Malignant
  # Otherwise use the other classifications
  return(obj)
}

# ---- 6. Neftel Module Scoring ----
score_neftel <- function(obj) {
  for (mod in names(NEFTEL_MODULES)) {
    genes <- intersect(NEFTEL_MODULES[[mod]], rownames(obj))
    if (length(genes) > 0) {
      obj <- AddModuleScore(obj, features = list(genes), name = paste0("Neftel_", mod))
      colnames(obj@meta.data)[ncol(obj@meta.data)] <- paste0("Neftel_", mod)
    }
  }
  return(obj)
}

# ---- 7. Compute WEBS ----
source("tools/WEBS/webs.R")

compute_webs_pipeline <- function(obj) {
  # UCell-based WEBS computation
  obj <- webs_compute(obj, 
                      writer_genes = WRITER_GENES,
                      eraser_genes = ERASER_GENES,
                      assay = "RNA", slot = "data")
  
  # Classify into LLS/LDTS/LHS by IQR
  obj <- webs_classify(obj, method = "iqr")
  
  return(obj)
}

# ---- 8. Region Type Classification ----
classify_region <- function(obj) {
  # Core score = mean(MKI67, TOP2A, PCNA)
  # Invasion score = mean(FN1, ITGA5, SPP1, MMP2)
  # Normal score = mean(RBFOX3, MAP2, GFAP)
  
  core_genes <- c("MKI67", "TOP2A", "PCNA")
  invasion_genes <- c("FN1", "ITGA5", "SPP1", "MMP2")
  normal_genes <- c("RBFOX3", "MAP2", "GFAP")
  
  get_mean_expr <- function(genes) {
    genes <- intersect(genes, rownames(obj))
    if (length(genes) == 0) return(rep(0, ncol(obj)))
    Matrix::colMeans(GetAssayData(obj, layer = "data")[genes, , drop = FALSE])
  }
  
  obj$Core_score <- get_mean_expr(core_genes)
  obj$Invasion_score <- get_mean_expr(invasion_genes)
  obj$Normal_score <- get_mean_expr(normal_genes)
  
  scores <- data.frame(
    Core = obj$Core_score,
    Invasion = obj$Invasion_score,
    Normal = obj$Normal_score
  )
  
  obj$Region_Type <- colnames(scores)[max.col(scores)]
  obj$Region_Type <- factor(obj$Region_Type, levels = c("Core", "Invasion", "Normal", "Intermediate"))
  
  return(obj)
}

# ---- 9. Full Pipeline ----
run_full_pipeline <- function(obj) {
  message("Step 1: QC...")
  obj <- run_qc(obj)
  
  message("Step 2: Preprocessing...")
  obj <- preprocess(obj)
  
  message("Step 3: Dimensionality reduction...")
  obj <- run_dimred(obj)
  
  message("Step 4: Cell type annotation...")
  obj <- annotate_cells(obj)
  
  message("Step 5: Neftel module scoring...")
  obj <- score_neftel(obj)
  
  message("Step 6: Compute WEBS...")
  obj <- compute_webs_pipeline(obj)
  
  message("Step 7: Region type classification...")
  obj <- classify_region(obj)
  
  message("Pipeline complete! Object has ", ncol(obj), " cells.")
  return(obj)
}

# ---- SAVE ----
# saveRDS(obj, "data/processed/GSE138794_webs_annotated.rds")
