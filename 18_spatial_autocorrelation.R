# =============================================================================
# 18. Spatial Autocorrelation of WEBS (Moran's I / Geary's c, permutation test)
# =============================================================================
# Tests whether high-WEBS spots cluster spatially in Visium slices beyond
# chance. Runs per slice on the spatial WEBS objects from 04_spatial.R.
# Input:  data/processed/spatial_<slice>.rds  (Seurat objects with $WEBS)
# Output: results/spatial_autocorrelation.csv
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(Seurat)
  library(spdep)
  library(dplyr)
})
set.seed(42)

N_PERM <- 1000

spatial_files <- list.files("data/processed", pattern = "^spatial_.*\\.rds$",
                            full.names = TRUE)
if (length(spatial_files) == 0) {
  stop("No spatial objects found in data/processed/. Run the 04_spatial.R driver first.")
}

run_autocorr <- function(obj, slice_id) {
  # Spot coordinates from the tissue image
  coords <- GetTissueCoordinates(obj)
  webs <- obj$WEBS
  # align
  common <- intersect(rownames(coords), names(webs))
  coords <- coords[common, ]
  webs <- webs[common]
  coords_mat <- as.matrix(coords[, 1:2])

  # k-nearest-neighbour weights (6 nearest, like Visium hex neighbours)
  knn <- knearneigh(coords_mat, k = 6)
  nb <- knn2nb(knn)
  lw <- nb2listw(nb, style = "W")

  moran <- moran.test(webs, lw)
  geary <- geary.test(webs, lw)

  # permutation p-value for Moran's I (label shuffling)
  I_obs <- moran$estimate["Moran I statistic"]
  I_perm <- replicate(N_PERM, {
    moran.test(sample(webs), lw)$estimate["Moran I statistic"]
  })
  p_perm <- (sum(I_perm >= I_obs) + 1) / (N_PERM + 1)

  data.frame(
    slice = slice_id,
    n_spots = length(webs),
    moran_I = I_obs,
    moran_p_analytic = moran$p.value,
    moran_p_permutation = p_perm,
    geary_c = geary$estimate["Geary C statistic"],
    geary_p = geary$p.value
  )
}

results <- lapply(spatial_files, function(f) {
  sid <- sub("^spatial_(.*)\\.rds$", "\\1", basename(f))
  message("Processing ", sid)
  obj <- readRDS(f)
  run_autocorr(obj, sid)
})
res_df <- do.call(rbind, results)

dir.create("results", showWarnings = FALSE)
write.csv(res_df, "results/spatial_autocorrelation.csv", row.names = FALSE)
print(res_df)

sink("results/spatial_autocorrelation_summary.txt")
cat("=== Spatial autocorrelation of WEBS across Visium slices ===\n\n")
print(res_df)
cat("\nInterpretation: Moran's I > 0 with permutation p < 0.05 indicates\n")
cat("significant spatial clustering of WEBS beyond random distribution.\n")
sink()

message("Done. See results/spatial_autocorrelation.csv")
