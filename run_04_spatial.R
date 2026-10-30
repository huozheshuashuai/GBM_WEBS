# =============================================================================
# Driver: Spatial transcriptomics on 6 GBM Visium slices (paper Figure 7)
# =============================================================================
# Input: data/raw/GSE194329/GSM*_spaceranger_out/ (extracted spaceranger dirs)
# Output: data/processed/spatial_<slice>.rds, results/spatial_summary.txt,
#         figures/fig_spatial_*.pdf
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(Seurat); library(ggplot2); library(dplyr); library(UCell) })
source("tools/WEBS/webs.R")

WRITERS <- c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP")
ERASERS <- c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")

slices <- list.dirs("data/raw/GSE194329", recursive = FALSE)
slices <- slices[grepl("GBM", slices)]

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("data/processed", showWarnings = FALSE)

results <- list()
for (sl in slices) {
  sid <- sub("_spaceranger_out$", "", basename(sl))
  sid <- sub("^GSM[0-9]+_", "", sid)
  out <- paste0("data/processed/spatial_", sid, ".rds")
  message("=== ", sid, " ===")

  inner <- list.dirs(sl, recursive = FALSE)
  data_dir <- inner[grepl("spaceranger_out$|GBM", basename(inner))][1]
  if (is.na(data_dir)) data_dir <- sl
  if (!file.exists(file.path(data_dir, "filtered_feature_bc_matrix.h5"))) {
    # search one level deeper
    cand <- list.files(data_dir, pattern = "filtered_feature_bc_matrix.h5",
                       recursive = TRUE, full.names = TRUE)
    if (length(cand) == 0) { message("no h5 found, skip"); next }
    data_dir <- dirname(cand[1])
  }

  sobj <- Load10X_Spatial(data.dir = data_dir,
                          filename = "filtered_feature_bc_matrix.h5",
                          assay = "Spatial", slice = sid, filter.matrix = TRUE)
  sobj <- SCTransform(sobj, assay = "Spatial", verbose = FALSE)

  # WEBS on SCT-normalized data
  sobj <- webs_compute(sobj, assay = "SCT", slot = "data")

  # Lactylation-proxy and metabolism scores (paper Methods)
  sct_data <- GetAssayData(sobj, assay = "SCT", layer = "data")
  metab_genes <- intersect(c("LDHA", "LDHB", "PKM", "SLC16A3", "ENO1", "GAPDH", "HK2", "PGK1"),
                           rownames(sct_data))
  sobj$Metabolism_score <- Matrix::colMeans(sct_data[metab_genes, , drop = FALSE])

  # dual-hotspot: lactylation-proxy (WEBS writer score) & metabolism > 75th pct
  w75 <- quantile(sobj$Writer_UCell, 0.75); m75 <- quantile(sobj$Metabolism_score, 0.75)
  sobj$dual_hotspot <- sobj$Writer_UCell > w75 & sobj$Metabolism_score > m75

  saveRDS(sobj, out)

  results[[sid]] <- data.frame(
    slice = sid, n_spots = ncol(sobj),
    mean_WEBS = mean(sobj$WEBS),
    n_dual_hotspot = sum(sobj$dual_hotspot)
  )

  p1 <- SpatialFeaturePlot(sobj, features = "WEBS", pt.size.factor = 1.6) +
    ggtitle(paste(sid, "WEBS"))
  pdf(paste0("figures/fig_spatial_", sid, ".pdf"), width = 6, height = 6)
  print(p1)
  dev.off()
}

res_df <- do.call(rbind, results)
write.csv(res_df, "results/spatial_summary.csv", row.names = FALSE)
print(res_df)
message("Done.")
