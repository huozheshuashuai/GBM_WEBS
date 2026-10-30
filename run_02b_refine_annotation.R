# =============================================================================
# Driver: refine cell-type annotation with original study labels
# =============================================================================
# Maps GSE138794_scRNA_Seq_cell_types.txt labels onto our merged object,
# defines malignant (Neoplastic) cells, recomputes WEBS stats on pure subsets.
# Input: data/processed/GSE138794_webs_annotated.rds
# Output: data/processed/GSE138794_final.rds, results/scrna_refined_summary.txt
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(Seurat); library(dplyr); library(ggplot2) })

obj <- readRDS("data/processed/GSE138794_webs_annotated.rds")

# ---- 1. Map original annotations ----
anno <- read.table(gzfile("data/raw/GSE138794/GSE138794_scRNA_Seq_cell_types.txt.gz"),
                   header = FALSE, sep = " ", col.names = c("cell_id", "original_type"),
                   fill = TRUE, stringsAsFactors = FALSE)
# our cell names: add.cell.ids used GSM sample names -> "GSMxxxx_Sample_barcode"
# annotation ids are "Sample_barcode" -> strip GSM prefix from ours
our_ids <- colnames(obj)
stripped <- sub("^GSM[0-9]+_", "", our_ids)
# barcodes may have trailing "-1"; annotation file may not. Harmonize both ways
anno$stripped_bc <- sub("-[0-9]+$", "", anno$cell_id)
our_bc <- sub("-[0-9]+$", "", stripped)
m <- match(our_bc, anno$stripped_bc)
obj$original_type <- anno$original_type[m]
cat("Annotated cells:", sum(!is.na(obj$original_type)), "/", ncol(obj), "\n")
print(table(obj$original_type, useNA = "ifany"))

# ---- 2. Malignant subset (Neoplastic) ----
obj$is_malignant <- obj$original_type == "Neoplastic_cell"
mal <- subset(obj, subset = is_malignant)
cat("Malignant (Neoplastic) cells:", ncol(mal), "\n")

# WEBS group by original cell type
sink("results/scrna_refined_summary.txt")
cat("=== Refined annotation summary ===\n\n")
cat("Total cells:", ncol(obj), "| annotated:", sum(!is.na(obj$original_type)),
    "| malignant:", ncol(mal), "\n\n")
cat("WEBS by original cell type (mean +/- sd):\n")
print(aggregate(WEBS ~ original_type, data = obj@meta.data,
                function(x) c(mean = mean(x), sd = sd(x))))
cat("\nWEBS group distribution among malignant cells:\n")
print(table(mal$WEBS_group))
cat("\nMalignant vs non-malignant WEBS (Wilcoxon):\n")
print(wilcox.test(mal$WEBS, obj$WEBS[!obj$is_malignant]))
cat("\nMalignant WEBS by sample:\n")
print(tapply(mal$WEBS, mal$sample, median))
sink()
cat(readLines("results/scrna_refined_summary.txt"), sep = "\n")

# ---- 3. Save ----
saveRDS(obj, "data/processed/GSE138794_final.rds")
saveRDS(mal, "data/processed/GSE138794_malignant.rds")
message("Saved GSE138794_final.rds and GSE138794_malignant.rds")
