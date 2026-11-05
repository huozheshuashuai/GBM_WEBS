# =============================================================================
# 65. GLASS recurrence: first-pair-per-patient sensitivity analysis
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr) })

lp <- read.csv("results/glass_recurrence_interval.csv")
sink("results/glass_firstpair_sensitivity.txt", split = TRUE)
cat("=== GLASS recurrence: first-pair-per-patient sensitivity ===\n\n")
cat("Full: n =", nrow(lp), "pairs,", length(unique(lp$case_barcode)), "patients\n")
ct <- cor.test(lp$WEBS_TP, lp$surgical_interval_mo, method = "spearman")
cat(sprintf("Full:        rho = %.3f, p = %.3g\n", ct$estimate, ct$p.value))

# keep one pair per patient (first recurrence only)
fp <- lp[!duplicated(lp$case_barcode), ]
ct2 <- cor.test(fp$WEBS_TP, fp$surgical_interval_mo, method = "spearman")
cat(sprintf("First-pair:  n = %d patients, rho = %.3f, p = %.3g\n",
            nrow(fp), ct2$estimate, ct2$p.value))
wt <- wilcox.test(fp$surgical_interval_mo ~ (fp$WEBS_TP > median(fp$WEBS_TP)))
cat(sprintf("First-pair median interval: low-WEBS %.0f mo vs high-WEBS %.0f mo, Wilcoxon p = %.3g\n",
            median(fp$surgical_interval_mo[fp$WEBS_TP <= median(fp$WEBS_TP)]),
            median(fp$surgical_interval_mo[fp$WEBS_TP > median(fp$WEBS_TP)]),
            wt$p.value))
sink()
cat(readLines("results/glass_firstpair_sensitivity.txt"), sep = "\n")
