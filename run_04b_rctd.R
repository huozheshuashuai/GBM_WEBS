# =============================================================================
# Driver: RCTD spatial deconvolution (spacexr) for 6 GBM Visium slices
# =============================================================================
# Reference: scRNA GSE138794 with original cell-type labels
# Input:  data/processed/GSE138794_final.rds, data/processed/spatial_*.rds
# Output: results/rctd_summary.txt, data/processed/spatial_<id>_rctd.rds
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(Seurat); library(spacexr); library(dplyr) })
set.seed(42)

# Reference: use annotated cells, top variable genes only (dense coercion safety)
ref <- readRDS("data/processed/GSE138794_final.rds")
ref <- subset(ref, cells = colnames(ref)[!is.na(ref$original_type)])
ref <- FindVariableFeatures(ref, nfeatures = 5000, verbose = FALSE)
hvgs <- VariableFeatures(ref)
ref_counts <- GetAssayData(ref, layer = "counts")[hvgs, ]
cell_types <- factor(ref$original_type)
names(cell_types) <- colnames(ref)

rctd_ref <- Reference(ref_counts, cell_types, n_max_cells = 6000)

spatial_files <- list.files("data/processed", pattern = "^spatial_.*\\.rds$",
                            full.names = TRUE)
spatial_files <- spatial_files[!grepl("rctd", spatial_files)]

sink("results/rctd_summary.txt")
cat("=== RCTD spatial deconvolution ===\n\n")
for (f in spatial_files) {
  sid <- sub("^spatial_(.*)\\.rds$", "\\1", basename(f))
  message("RCTD: ", sid)
  sobj <- readRDS(f)
  coords_df <- GetTissueCoordinates(sobj)
  sp_counts <- GetAssayData(sobj, assay = "Spatial", layer = "counts")
  nUMI <- colSums(sp_counts)
  puck <- SpatialRNA(coords = coords_df[, 1:2], counts = sp_counts, nUMI = nUMI)
  my_rctd <- create.RCTD(puck, rctd_ref, max_cores = 4)
  my_rctd <- run.RCTD(my_rctd, doublet_mode = "doublet")

  w <- my_rctd@results$weights
  sobj$rctd_malignant <- w[colnames(sobj), "Neoplastic_cell"]
  sobj$rctd_myeloid <- w[colnames(sobj), "Myeloid"]

  ct <- cor.test(sobj$rctd_malignant, sobj$WEBS, method = "spearman")
  ct2 <- cor.test(sobj$rctd_myeloid, sobj$WEBS, method = "spearman")
  cat(sid, ": WEBS vs malignant fraction rho =", round(ct$estimate, 3),
      " p =", signif(ct$p.value, 3),
      " | vs myeloid rho =", round(ct2$estimate, 3),
      " p =", signif(ct2$p.value, 3), "\n")
  saveRDS(sobj, sub("\\.rds$", "_rctd.rds", f))
}
sink()
cat(readLines("results/rctd_summary.txt"), sep = "\n")
message("Done.")
