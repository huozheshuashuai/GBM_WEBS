# =============================================================================
# 52. Fix panel issues + assemble composite Figures 13 & 14
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork); library(survival); library(survminer) })

PURPLE <- "#7B4EA3"; TEAL <- "#3D9B8C"; RED <- "#d62728"; BLUE <- "#1f77b4"

# ---------- panel 13a: splicing volcano ----------
res <- read.csv("results/splicing_webs_corr.csv")
res$sig_flag <- ifelse(res$q < 0.05 & abs(res$rho) > 0.2, "Sig", "NS")
nsig <- sum(res$sig_flag == "Sig")
p13a <- ggplot(res, aes(rho, -log10(p), color = sig_flag)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Sig" = PURPLE, "NS" = "grey80")) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey40") +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = "Spearman rho (PSI vs WEBS)", y = "-log10(p)",
       title = paste0("AS events vs WEBS (", nsig, " sig)"))

# ---------- panel 13b: event types ----------
sig <- res[res$sig_flag == "Sig", ]
tt <- as.data.frame(table(sig$splice_type)); colnames(tt) <- c("type", "count")
p13b <- ggplot(tt, aes(reorder(type, -count), count, fill = type)) +
  geom_bar(stat = "identity") + scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = "Event type", y = "Significant events", title = "Event-type composition")

# ---------- panel 13c: enrichment (fixed title) ----------
ego <- read.csv("results/splicing_enrichment.csv")
ego <- ego[order(ego$p.adjust), ][1:8, ]
ego$Description <- factor(ego$Description, levels = rev(ego$Description))
ego$neglog <- -log10(ego$p.adjust)
p13c <- ggplot(ego, aes(neglog, Description, size = Count, color = p.adjust)) +
  geom_point() + scale_color_gradient(low = PURPLE, high = "grey60") +
  theme_bw(base_size = 9) +
  theme(legend.position = "right", legend.key.size = unit(0.35, "cm"),
        axis.text.y = element_text(size = 7.5)) +
  labs(x = "-log10(adj. p)", y = NULL, size = "Genes", color = "adj. p",
       title = "GO enrichment of AS genes")

# ---------- panel 13d: phospho volcano (spliceosome on TOP) ----------
ph <- read.csv("results/phospho_limma_LHS_vs_LLS.csv")
ph$cat <- ifelse(ph$adj.P.Val < 0.05 & ph$is_spliceosome, "Spliceosome sig",
          ifelse(ph$adj.P.Val < 0.05, "Other sig",
          ifelse(ph$is_spliceosome, "Spliceosome NS", "NS")))
ph <- ph[order(ph$cat != "Spliceosome sig"), ]   # spliceosome sig drawn last
ph$pt_size <- ifelse(ph$cat == "Spliceosome sig", 1.6, 0.5)
p13d <- ggplot(ph, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = cat, size = pt_size), alpha = 0.55) +
  scale_size_identity() +
  scale_color_manual(values = c("Spliceosome sig" = PURPLE, "Other sig" = "grey80",
                                "Spliceosome NS" = TEAL, "NS" = "grey92")) +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = "logFC (LHS vs LLS)", y = "-log10(p)",
       title = "Phosphosites: 198/319 spliceosome sig")

# ---------- panel 13e: kinase correlations ----------
webs <- read.csv("results/cptac_webs.csv")
prot <- read.delim("E:/KIMI/data/HS_CPTAC_GBM_proteome_mssm_per_gene.cct", check.names = FALSE)
pm <- as.matrix(sapply(prot[, -1], as.numeric)); rownames(pm) <- prot[[1]]
kinases <- c("SRPK1","SRPK2","CLK1","CLK2","CLK3","CLK4","CDK7","CDK9","CDK12","CDK13")
kinases <- intersect(kinases, rownames(pm))
kk <- intersect(colnames(pm), webs$sample)
kdf <- data.frame(kinase = kinases, rho = NA, p = NA)
for (k in kinases) {
  ct <- suppressWarnings(cor.test(pm[k, kk], webs$WEBS[match(kk, webs$sample)], method = "spearman"))
  kdf$rho[kdf$kinase == k] <- ct$estimate; kdf$p[kdf$kinase == k] <- ct$p.value
}
kdf$sig <- ifelse(kdf$p < 0.05, "p<0.05", "n.s.")
p13e <- ggplot(kdf, aes(reorder(kinase, rho), rho, fill = sig)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("p<0.05" = PURPLE, "n.s." = "grey70")) +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = NULL, y = "Spearman rho vs WEBS", title = "Splicing kinases vs WEBS")

# ---------- panel 13f: virtual KO (absolute sensitivity ratio) ----------
ko <- read.csv("results/virtual_ko_abs_summary.csv")
ko$type <- ifelse(ko$ko %in% c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP"), "Writer", "Eraser")
p13f <- ggplot(ko, aes(reorder(ko, ratio), ratio, fill = type)) +
  geom_bar(stat = "identity") + coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("Writer" = PURPLE, "Eraser" = TEAL)) +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = NULL, y = "Relative |shift| (spliceosome / other)",
       title = "Virtual KO: spliceosome sensitivity")

# ---------- assemble Figure 13 ----------
fig13 <- (p13a | p13b | p13c) / (p13d | p13e | p13f) +
  plot_annotation(tag_levels = "a")
ggsave("figures_final/Figure_13.pdf", fig13, width = 14, height = 9)
ggsave("figures_final/Figure_13.png", fig13, width = 14, height = 9, dpi = 200)

# ---------- panel 14a: GLASS KM ----------
prim <- read.csv("results/glass_primary_scores.csv")
prim$risk_group <- factor(ifelse(prim$risk > median(prim$risk), "High", "Low"), levels = c("Low", "High"))
fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = prim)
p14a <- ggsurvplot(fit, data = prim, pval = TRUE, legend.labs = c("Low risk", "High risk"),
                   palette = c(BLUE, RED), xlab = "Time (months)", ylab = "Overall survival",
                   title = "GLASS primary (n=107)", legend = "top",
                   font.main = 11, font.x = 10, font.y = 10)$plot

# ---------- panel 14b: recurrence time ----------
lp <- read.csv("results/glass_recurrence_interval.csv")
ct <- cor.test(lp$WEBS_TP, lp$surgical_interval_mo, method = "spearman")
p14b <- ggplot(lp, aes(WEBS_TP, surgical_interval_mo)) +
  geom_point(alpha = 0.5, color = PURPLE, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  theme_bw(base_size = 10) +
  labs(x = "Primary tumor WEBS", y = "Months to recurrence",
       title = paste0("GLASS recurrence (rho=", round(ct$estimate, 2),
                      ", p=", signif(ct$p.value, 1), ")"))

# ---------- panel 14c/d: pan-cancer (GBM highlighted) ----------
pc <- read.csv("results/pancancer_33_summary.csv")
pc$cat <- ifelse(pc$cancer == "GBM", "GBM", ifelse(pc$q_webs_splice < 0.05, "q<0.05", "n.s."))
pc$cancer <- factor(pc$cancer, levels = pc$cancer[order(pc$rho_webs_splice)])
p14c <- ggplot(pc, aes(rho_webs_splice, cancer, color = cat, size = cat)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("GBM" = RED, "q<0.05" = PURPLE, "n.s." = "grey65")) +
  scale_size_manual(values = c("GBM" = 3.5, "q<0.05" = 2, "n.s." = 2)) +
  theme_bw(base_size = 9) + theme(legend.position = "none") +
  labs(x = "Spearman rho", y = NULL, title = "WEBS-spliceosome (33 TCGA)")

pc2 <- pc[!is.na(pc$sig_HR), ]
pc2$cat <- ifelse(pc2$cancer == "GBM", "GBM", ifelse(pc2$sig_p < 0.05, "p<0.05", "n.s."))
pc2$cancer <- factor(pc2$cancer, levels = pc2$cancer[order(pc2$sig_HR)])
p14d <- ggplot(pc2, aes(sig_HR, cancer, color = cat, size = cat)) +
  geom_point() +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c("GBM" = RED, "p<0.05" = PURPLE, "n.s." = "grey65")) +
  scale_size_manual(values = c("GBM" = 3.5, "p<0.05" = 2, "n.s." = 2)) +
  theme_bw(base_size = 9) + theme(legend.position = "none") +
  labs(x = "Hazard ratio", y = NULL, title = "Six-gene signature HR (33 TCGA)")

# ---------- assemble Figure 14 ----------
fig14 <- (p14a | p14b) / (p14c | p14d) +
  plot_annotation(tag_levels = "a")
ggsave("figures_final/Figure_14.pdf", fig14, width = 13, height = 11)
ggsave("figures_final/Figure_14.png", fig14, width = 13, height = 11, dpi = 200)
cat("Done.\n")
