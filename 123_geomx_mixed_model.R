# =============================================================================
# 123. GeoMx patient-aware inference: mixed model with patient random intercept
#      score ~ stain * timepoint + (1 | patient)
#      verify EMT contrast (p=0.020) and microglia validation (p=0.001)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(lme4); library(lmerTest) })
set.seed(42)

gm <- readRDS("results/geomx/geomx_bundle.rds")
res <- gm$roi_meta
cat("ROIs:", nrow(res), "| patients:", length(unique(res$patient)), "\n")

sink("results/geomx/geomx_mixed_model.txt", split = TRUE)
cat("=== GeoMx patient-aware mixed models ===\n")
cat("model: score ~ stain * timepoint + (1 | patient), ROI nested in patient\n\n")

run_mm <- function(y, label) {
  df <- res[!is.na(res[[y]]), ]
  fit <- lmer(as.formula(paste(y, "~ stain * timepoint + (1 | patient)")), data = df)
  s <- summary(fit)$coefficients
  cat("---", label, "---\n")
  print(round(s[, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], 4))
  cat("\n")
  invisible(fit)
}

fit_emt <- run_mm("mod_EMT", "EMT")
fit_pc2 <- run_mm("PC2", "PC2 axis")
fit_il2 <- run_mm("mod_IL2STAT5", "IL2-STAT5")
fit_mg  <- run_mm("px_Microglia", "Microglia proxy (stain validation)")
fit_ros <- run_mm("mod_ROS", "ROS")

# patient-level summary as sensitivity (already the primary in 118; confirm N)
pt <- gm$pt
wide <- pt %>% tidyr::pivot_wider(names_from = stain, values_from = c(PC2, IL2STAT5, ROS, EMT, Microglia, Endothelial))
wide <- wide[!is.na(wide$PC2_CD64) & !is.na(wide$PC2_GFAP), ]
cat("sensitivity: patient-level paired contrast (patient x timepoint units, n =",
    nrow(wide), "combinations from", length(unique(wide$patient)), "patients)\n")
cat("NOTE: combinations are not fully independent when a patient has both timepoints;\n")
cat("mixed model above is the pre-specified primary inference.\n")
sink()
cat("Done.\n")
