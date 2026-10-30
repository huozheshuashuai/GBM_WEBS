# =============================================================================
# Driver: GSE138794 scRNA-seq full pipeline (QC -> WEBS -> Region types)
# =============================================================================
# Usage: Rscript run_02_scrna.R
# Working directory: GBM_WEBS_Project/
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
set.seed(42)

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(UCell)
})

source("tools/WEBS/webs.R")

WRITER_GENES <- c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP")
ERASER_GENES <- c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")

# ---- 1. Locate extracted per-sample files and reorganize into 10x dirs ----
raw_dir <- "data/raw/GSE138794/extracted"
prep_dir <- "data/raw/GSE138794/tenx"
dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)

all_files <- list.files(raw_dir, recursive = FALSE)
# scRNA samples: have features.tsv.gz or genes.tsv.gz (exclude snATAC which use peaks.bed.gz)
feature_files <- grep("features\\.tsv\\.gz$|genes\\.tsv\\.gz$", all_files, value = TRUE)
# exclude any snATAC-labelled
feature_files <- grep("snATAC", feature_files, value = TRUE, invert = TRUE)

message("Found ", length(feature_files), " scRNA samples")

samples <- sub("_(features|genes)\\.tsv\\.gz$", "", feature_files)

for (s in samples) {
  sdir <- file.path(prep_dir, s)
  dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
  ff <- list.files(raw_dir, pattern = paste0("^", s, "_"), full.names = TRUE)
  for (f in ff) {
    base <- basename(f)
    newname <- sub(paste0("^", s, "_"), "", base)
    # normalize genes.tsv.gz -> features.tsv.gz
    if (newname == "genes.tsv.gz") newname <- "features.tsv.gz"
    target <- file.path(sdir, newname)
    if (!file.exists(target)) file.copy(f, target)
  }
}

# ---- 2. Read and merge all samples ----
objs <- list()
for (s in samples) {
  sdir <- file.path(prep_dir, s)
  needed <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  if (!all(file.exists(file.path(sdir, needed)))) {
    message("Skipping ", s, " (incomplete 10x files)")
    next
  }
  # features.tsv may have 1 column (older samples) or 2-3 (standard 10x)
  first_line <- readLines(gzfile(file.path(sdir, "features.tsv.gz")), n = 1)
  n_cols <- length(strsplit(first_line, "\t")[[1]])
  counts <- Read10X(data.dir = sdir, gene.column = min(2, n_cols))
  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  sample_name <- sub("^GSM[0-9]+_", "", s)
  obj <- CreateSeuratObject(counts = counts, project = sample_name,
                            min.cells = 3, min.features = 200)
  obj$sample <- sample_name
  objs[[s]] <- obj
  message("Loaded ", s, ": ", ncol(obj), " cells")
}

message("Merging ", length(objs), " samples...")
obj <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))
rm(objs); gc()

# Seurat v5 keeps per-sample layers split; join them for downstream ops
obj <- JoinLayers(obj)
message("Merged object: ", nrow(obj), " features x ", ncol(obj), " cells")
message("RBFOX3 present: ", "RBFOX3" %in% rownames(obj))

# ---- 3. QC ----
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

n_before <- ncol(obj)
obj <- subset(obj, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 20)
message("QC: ", n_before, " -> ", ncol(obj), " cells")

# ---- 4. Normalize, variable features, scaling ----
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)

# ---- 5. PCA / clustering / UMAP ----
obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.8)
obj <- RunUMAP(obj, dims = 1:30)

# ---- 6. Cell type annotation (marker-based max score) ----
CELL_TYPE_MARKERS <- list(
  Neuron = c("RBFOX3", "MAP2"),
  Astrocyte = c("GFAP"),
  Macrophage_Microglia = c("PTPRC", "CD14", "CSF1R", "C1QA", "AIF1"),
  Endothelial = c("PECAM1"),
  Pericyte = c("PDGFRB"),
  Oligodendrocyte = c("OLIG1", "OLIG2", "MBP"),
  Malignant = c("SOX2", "EGFR", "OLIG2")
)

for (ct in names(CELL_TYPE_MARKERS)) {
  markers <- intersect(CELL_TYPE_MARKERS[[ct]], rownames(obj))
  if (length(markers) > 0) {
    obj <- AddModuleScore(obj, features = list(markers), name = paste0("CT_", ct))
  }
}
score_cols <- paste0("CT_", names(CELL_TYPE_MARKERS), "1")
score_cols <- intersect(score_cols, colnames(obj@meta.data))
obj$cell_type <- sub("^CT_(.*)1$", "\\1", score_cols[max.col(obj@meta.data[, score_cols])])

# ---- 7. Neftel module scoring (top markers per module) ----
NEFTEL <- list(
  MES = c("CD44", "VIM", "S100A11", "ANXA1", "LGALS1", "BCL2A1", "CD83", "CXCL8", "AREG", "CD74"),
  AC  = c("APOE", "GFAP", "S100B", "HEPN1", "GJA1", "SPARCL1", "CLU", "AGT", "FABP7", "SLC1A3"),
  NPC = c("SOX4", "STMN1", "STMN2", "BEX1", "SOX11", "DLL3", "MARCKSL1", "GAP43", "INA", "NNAT"),
  OPC = c("SOX10", "PLP1", "BCAS1", "GPR17", "PTPRZ1", "CA2", "MYRF", "MAG", "MOG", "APOD")
)
for (mod in names(NEFTEL)) {
  genes <- intersect(NEFTEL[[mod]], rownames(obj))
  obj <- AddModuleScore(obj, features = list(genes), name = paste0("Neftel_", mod))
  colnames(obj@meta.data)[ncol(obj@meta.data)] <- paste0("Neftel_", mod)
}

# ---- 8. WEBS ----
obj <- webs_compute(obj, writer_genes = WRITER_GENES, eraser_genes = ERASER_GENES)
obj <- webs_classify(obj, method = "iqr")

# ---- 9. Region type classification ----
get_mean_expr <- function(o, genes) {
  genes <- intersect(genes, rownames(o))
  if (length(genes) == 0) return(rep(0, ncol(o)))
  Matrix::colMeans(GetAssayData(o, layer = "data")[genes, , drop = FALSE])
}
obj$Core_score <- get_mean_expr(obj, c("MKI67", "TOP2A", "PCNA"))
obj$Invasion_score <- get_mean_expr(obj, c("FN1", "ITGA5", "SPP1", "MMP2"))
obj$Normal_score <- get_mean_expr(obj, c("RBFOX3", "MAP2", "GFAP"))
scores <- data.frame(Core = obj$Core_score, Invasion = obj$Invasion_score, Normal = obj$Normal_score)
obj$Region_Type <- colnames(scores)[max.col(scores)]

# Map paper nomenclature: Invasion -> Margin-like
obj$Region_Type <- factor(obj$Region_Type, levels = c("Core", "Invasion", "Normal"),
                          labels = c("Core-like", "Margin-like", "Normal-like"))

# ---- 10. Hypoxia scores ----
hallmark_hypoxia <- c("ACOT7","ADM","ALDOA","AMPD3","ATF3","BHLHE40","BNIP3","BNIP3L",
  "CA9","CCNG2","CDKN1A","CDKN1B","CDKN1C","CHCHD2","CITED2","CP","CXCR4","DDIT4","ENO1",
  "ERRFI1","FAM162A","GPI","HK2","HMOX1","HSPA5","IGFBP3","LDHA","MIF","NDRG1","P4HA1",
  "PFKFB3","PFKP","PGK1","PKM","PLAUR","SLC2A1","STC2","TPI1","VEGFA","ZFP36")
hallmark_hypoxia <- intersect(hallmark_hypoxia, rownames(obj))
obj <- AddModuleScore(obj, features = list(hallmark_hypoxia), name = "Hallmark_Hypoxia")
colnames(obj@meta.data)[ncol(obj@meta.data)] <- "Hallmark_Hypoxia"
obj$CA9_VEGFA <- get_mean_expr(obj, c("CA9", "VEGFA"))

# ---- 11. Save ----
dir.create("data/processed", showWarnings = FALSE)
saveRDS(obj, "data/processed/GSE138794_webs_annotated.rds")

# ---- 12. Summary report ----
sink("results/scrna_pipeline_summary.txt")
cat("=== GSE138794 scRNA-seq pipeline summary ===\n\n")
cat("Total cells after QC:", ncol(obj), "\n")
cat("Samples:", length(unique(obj$sample)), "\n\n")
cat("Cell type distribution:\n"); print(table(obj$cell_type))
cat("\nWEBS group distribution:\n"); print(table(obj$WEBS_group))
cat("\nRegion type distribution:\n"); print(table(obj$Region_Type))
cat("\nWEBS by cell type (median):\n")
print(tapply(obj$WEBS, obj$cell_type, median))
cat("\nWEBS by region type (median):\n")
print(tapply(obj$WEBS, obj$Region_Type, median))
cat("\nMargin-like vs Core-like WEBS (Wilcoxon):\n")
m <- obj$WEBS[obj$Region_Type == "Margin-like"]; c2 <- obj$WEBS[obj$Region_Type == "Core-like"]
print(wilcox.test(m, c2))
cat("\nSpearman WEBS vs Hallmark_Hypoxia:\n")
print(cor.test(obj$WEBS, obj$Hallmark_Hypoxia, method = "spearman"))
cat("\nSpearman WEBS vs CA9/VEGFA:\n")
print(cor.test(obj$WEBS, obj$CA9_VEGFA, method = "spearman"))
sink()

message("Done. Object saved to data/processed/GSE138794_webs_annotated.rds")
