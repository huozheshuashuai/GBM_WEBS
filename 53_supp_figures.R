# =============================================================================
# 53. New supplementary figures S11-S19 (v2 submission package)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(ggplot2); library(patchwork) })

OUT <- "submission/Supplementary/figures/pdf"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# --- S12: treatment interaction (-log10 p barplot, from 38_confounders results) ---
ti <- data.frame(
  cohort = rep(c("CGGA-693", "CGGA-325"), each = 2),
  treatment = rep(c("Temozolomide", "Radiotherapy"), 2),
  p = c(0.00239, 0.000454, 0.000772, 0.0258))
ti$neglog <- -log10(ti$p)
p_ti <- ggplot(ti, aes(treatment, neglog, fill = cohort)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("CGGA-693" = "#7B4EA3", "CGGA-325" = "#3D9B8C")) +
  theme_bw(base_size = 12) +
  labs(title = "Risk x treatment interaction (Cox models)",
       x = NULL, y = "-log10(interaction p)", fill = NULL)
ggsave(file.path(OUT, "figS_treatment_interaction.pdf"), p_ti, width = 6, height = 4.5)

cat("Done R part.\n")
