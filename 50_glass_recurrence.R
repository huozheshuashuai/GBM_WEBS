# =============================================================================
# 50. GLASS: does primary-tumor WEBS predict time to recurrence?
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(survival); library(ggplot2); library(dplyr) })

DATA <- "E:/KIMI/data/GLASS_extracted"

tpm <- read.delim(file.path(DATA, "gene_tpm_matrix_all_samples.tsv"), check.names = FALSE)
expr <- as.matrix(tpm[, -1]); rownames(expr) <- tpm[[1]]
expr <- log2(expr + 1)

WRITERS <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP")
ERASERS <- c("HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")
webs <- colMeans(expr[intersect(WRITERS, rownames(expr)), , drop = FALSE]) -
        colMeans(expr[intersect(ERASERS, rownames(expr)), , drop = FALSE])
webs_df <- data.frame(aliquot = names(webs), WEBS = as.numeric(webs))

pairs <- read.csv(file.path(DATA, "analysis_rnaseq_pairs.csv"))
pairs$a_dot <- gsub("-", ".", pairs$tumor_barcode_a)
pairs$b_dot <- gsub("-", ".", pairs$tumor_barcode_b)
lp <- pairs[pairs$comparison_type == "longitudinal", ]
lp$WEBS_TP <- webs_df$WEBS[match(lp$a_dot, webs_df$aliquot)]
lp$WEBS_R1 <- webs_df$WEBS[match(lp$b_dot, webs_df$aliquot)]
lp$surgical_interval_mo <- as.numeric(lp$surgical_interval_mo)
lp <- lp[!is.na(lp$WEBS_TP) & !is.na(lp$surgical_interval_mo) & lp$surgical_interval_mo > 0, ]

sink("results/glass_recurrence_time.txt", split = TRUE)
cat("=== GLASS: primary WEBS vs time to recurrence ===\n")
cat("Pairs with interval:", nrow(lp), "\n")
ct <- cor.test(lp$WEBS_TP, lp$surgical_interval_mo, method = "spearman")
cat("Spearman rho =", round(ct$estimate, 3), " p =", signif(ct$p.value, 3), "\n")
# median split
lp$grp <- factor(ifelse(lp$WEBS_TP > median(lp$WEBS_TP), "High WEBS", "Low WEBS"),
                 levels = c("Low WEBS", "High WEBS"))
wt <- wilcox.test(surgical_interval_mo ~ grp, data = lp)
cat("Median interval Low:", median(lp$surgical_interval_mo[lp$grp == "Low WEBS"]),
    " High:", median(lp$surgical_interval_mo[lp$grp == "High WEBS"]),
    " Wilcoxon p =", signif(wt$p.value, 3), "\n")
# WEBS change vs interval: does longer interval associate with bigger WEBS jump?
lp$dWEBS <- lp$WEBS_R1 - lp$WEBS_TP
ok <- !is.na(lp$dWEBS)
ct2 <- cor.test(lp$surgical_interval_mo[ok], lp$dWEBS[ok], method = "spearman")
cat("Interval vs WEBS change: rho =", round(ct2$estimate, 3),
    " p =", signif(ct2$p.value, 3), "\n")
sink()

p <- ggplot(lp, aes(WEBS_TP, surgical_interval_mo)) +
  geom_point(alpha = 0.6, color = "#7B4EA3") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  theme_bw(base_size = 12) +
  labs(title = "GLASS: primary WEBS vs time to recurrence",
       subtitle = paste0("Spearman rho = ", round(ct$estimate, 2),
                         ", p = ", signif(ct$p.value, 2), " (n=", nrow(lp), ")"),
       x = "Primary tumor WEBS", y = "Surgical interval (months)")
ggsave("figures_final/fig_glass_recurrence_time.pdf", p, width = 6, height = 5)
write.csv(lp, "results/glass_recurrence_interval.csv", row.names = FALSE)
cat("Done.\n")
