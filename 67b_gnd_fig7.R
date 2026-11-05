# 67b. Finish GND Figure 7 (GLASS + pharmaco composite)
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork) })
PURPLE <- "#7B4EA3"; TEAL <- "#3D9B8C"

lp <- read.csv("results/glass_recurrence_interval.csv")
p7a <- ggplot(lp, aes(WEBS_TP, surgical_interval_mo)) +
  geom_point(alpha = 0.5, color = PURPLE, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", formula = y ~ x) +
  theme_bw(base_size = 10) +
  labs(x = "Primary tumor WEBS", y = "Months to recurrence",
       title = "GLASS: recurrence interval (rho=0.33)")

lo <- read.csv("results/glass_longitudinal_pairs.csv")
plot_df <- data.frame(case = rep(lo$case_barcode, 2),
                      stage = factor(rep(c("Primary", "Recurrence"), each = nrow(lo)),
                                     levels = c("Primary", "Recurrence")),
                      WEBS = c(lo$WEBS_TP, lo$WEBS_R1))
p7b <- ggplot(plot_df, aes(stage, WEBS, group = case)) +
  geom_line(alpha = 0.25, color = "grey50") +
  geom_point(aes(color = stage), size = 1.2) +
  scale_color_manual(values = c("Primary" = TEAL, "Recurrence" = PURPLE)) +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = NULL, y = "WEBS", title = "GLASS longitudinal pairs")

gd <- read.csv("results/gdsc2_drug_webs_corr.csv")
gd <- gd[!is.na(gd$q) & gd$q < 0.1, ]
top_drugs <- gd[order(gd$spearman_r), ][1:10, ]
top_drugs$DRUG_NAME <- factor(top_drugs$DRUG_NAME, levels = rev(top_drugs$DRUG_NAME))
p7c <- ggplot(top_drugs, aes(DRUG_NAME, spearman_r, fill = spearman_r < 0)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("TRUE" = PURPLE, "FALSE" = "grey70")) +
  theme_bw(base_size = 9) + theme(legend.position = "none") +
  labs(x = NULL, y = "Spearman rho (WEBS vs ln IC50)",
       title = "Top WEBS-sensitizing drugs (GDSC2, q<0.1)")

dd <- read.csv("results/depmap_webs_dependency.csv")
print(colnames(dd))
dd$gene <- dd[[1]]
for (cand in c("MED1","SF3B1","SRSF1","HDAC2","SIRT2")) {
  if (cand %in% dd$gene) next
}
key <- dd[dd$gene %in% c("MED1","SF3B1","SRSF1","HDAC2","SIRT2"), ]
rcol <- grep("rho|r\\b|cor", colnames(key), value = TRUE, ignore.case = TRUE)[1]
pcol <- grep("^q|adj", colnames(key), value = TRUE, ignore.case = TRUE)[1]
key$gene <- factor(key$gene, levels = key$gene[order(key[[rcol]])])
p7d <- ggplot(key, aes(gene, .data[[rcol]], fill = .data[[pcol]] < 0.05)) +
  geom_bar(stat = "identity") + coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("TRUE" = PURPLE, "FALSE" = "grey70")) +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = NULL, y = "Correlation with dependency (Chronos)",
       title = "DepMap dependency (glioma lines)")

fig7 <- (p7a | p7b) / (p7c | p7d) + plot_annotation(tag_levels = "a")
ggsave("figures_GND/Figure_7.pdf", fig7, width = 12, height = 9)
ggsave("figures_GND/Figure_7.png", fig7, width = 12, height = 9, dpi = 200)
cat("Figure 7 done.\n")
