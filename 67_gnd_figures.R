# =============================================================================
# 67. Genes & Diseases main figures (7F + 1T layout)
#     F3: CPTAC spliceosome accumulation (NEW - core result had no main figure)
#     F4: PTM layer (acetyl + phospho + kinases + virtual KO)
#     F5: Splicing output (+ LGG replication)
#     F6: Clinical validation (4-cohort KM)
#     F7: GLASS recurrence + pharmacogenomics
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork); library(survival); library(survminer) })
PURPLE <- "#7B4EA3"; TEAL <- "#3D9B8C"; RED <- "#d62728"; BLUE <- "#1f77b4"
dir.create("figures_GND", showWarnings = FALSE)

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
webs <- read.csv("results/cptac_webs.csv")

# ============ Figure 3: CPTAC spliceosome accumulation ============
# (a) volcano with spliceosome highlighted
lim <- read.csv("results/cptac_limma_results.csv")
if (!"is_spliceosome" %in% colnames(lim)) {
  gene_col <- colnames(lim)[sapply(lim, function(x) any(x %in% splice))][1]
  lim$is_spliceosome <- lim[[gene_col]] %in% splice
}
lim$cat <- ifelse(lim$is_spliceosome & lim$adj.P.Val < 0.05 & lim$logFC > 0, "Spliceosome up",
           ifelse(lim$is_spliceosome, "Spliceosome", "Other"))
lim <- lim[order(lim$cat != "Spliceosome up"), ]
lim$sz <- ifelse(lim$cat == "Spliceosome up", 2.0, 0.4)
lim$shp <- ifelse(lim$cat == "Spliceosome up", 17, 16)
p3a <- ggplot(lim, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = cat, size = sz, shape = shp), alpha = 0.65) +
  scale_shape_identity() +
  scale_size_identity() +
  scale_color_manual(values = c("Spliceosome up" = PURPLE, "Spliceosome" = TEAL, "Other" = "grey85")) +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = "logFC (high- vs low-WEBS)", y = "-log10(p)",
       title = "CPTAC proteome: 107/111 spliceosome up")

# (b) scatter protein WEBS vs spliceosome protein score
prot <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct", check.names = FALSE)
pm <- as.matrix(sapply(prot[, -1], as.numeric)); rownames(pm) <- prot[[1]]
cs <- intersect(colnames(pm), webs$sample)
sp <- intersect(splice, rownames(pm))
score <- colMeans(pm[sp, cs, drop = FALSE], na.rm = TRUE)
df3b <- data.frame(WEBS = webs$WEBS[match(cs, webs$sample)], spliceosome = score)
p3b <- ggplot(df3b, aes(WEBS, spliceosome)) +
  geom_point(color = PURPLE, alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey30", se = TRUE) +
  theme_bw(base_size = 10) +
  labs(x = "Protein-level WEBS", y = "Spliceosome protein score",
       title = "rho = 0.79 (95% CI 0.68-0.87)")

# (c) per-gene RNA vs protein rho
pg <- read.csv("results/uncoupling_per_gene.csv")
p3c <- ggplot(pg, aes(rho_mrna, rho_protein)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 1.8, alpha = 0.7, color = PURPLE) +
  theme_bw(base_size = 10) +
  labs(x = "mRNA-level rho", y = "Protein-level rho",
       title = "107/109 proteins above diagonal")

# (d) threshold robustness
td <- data.frame(cutoff = c("20%", "25%", "30%", "median"),
                 up = c(110, 110, 111, 109))
p3d <- ggplot(td, aes(cutoff, up)) +
  geom_bar(stat = "identity", fill = PURPLE) +
  geom_hline(yintercept = 111, linetype = "dashed", color = "grey40") +
  ylim(0, 115) +
  theme_bw(base_size = 10) +
  labs(x = "High/low-WEBS cutoff", y = "Spliceosome proteins up",
       title = "Threshold robustness (109-111/111)")

fig3 <- (p3a | p3b) / (p3c | p3d) + plot_annotation(tag_levels = "a")
ggsave("figures_GND/Figure_3.pdf", fig3, width = 10, height = 8)
ggsave("figures_GND/Figure_3.png", fig3, width = 10, height = 8, dpi = 200)

# ============ Figure 4: PTM layer ============
# (a) acetyl heatmap (reload from generated pdf is complex; rebuild quickly)
tt <- read.csv("results/acetyl_sites_limma.csv")
tt$gene <- sub("-.*$", "", tt$site)
top <- tt %>% filter(gene %in% splice, adj.P.Val < 0.05) %>% arrange(desc(abs(t))) %>% head(30)
ac <- read.delim("E:/KIMI/data/HS_CPTAC_GBM_acetylome_mssm_site.cct", check.names = FALSE)
sub <- ac[match(top$site, ac[[1]]), ]
am <- as.matrix(sapply(sub[, -1], as.numeric)); rownames(am) <- top$site
common <- intersect(colnames(am), webs$sample)
am <- am[, common]; w2 <- webs[match(common, webs$sample), ]
am <- am[, order(w2$WEBS)]
library(tidyr)
long <- pivot_longer(as.data.frame(am) %>% mutate(site = rownames(am)), -site,
                     names_to = "sample", values_to = "val")
long$site <- factor(long$site, levels = rev(top$site))
long$sample <- factor(long$sample, levels = colnames(am))
p4a <- ggplot(long, aes(sample, site, fill = val)) +
  geom_tile() +
  scale_fill_gradient2(low = TEAL, mid = "white", high = PURPLE, midpoint = 0, na.value = "grey90") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6)) +
  labs(x = "Samples ordered by WEBS", y = NULL, fill = "Acetyl.",
       title = "Spliceosome acetylation (OR=7.36; corrected 2.29)")

# (b) phospho volcano
ph <- read.csv("results/phospho_limma_LHS_vs_LLS.csv")
ph$cat <- ifelse(ph$adj.P.Val < 0.05 & ph$is_spliceosome, "Spliceosome sig",
          ifelse(ph$adj.P.Val < 0.05, "Other sig",
          ifelse(ph$is_spliceosome, "Spliceosome NS", "NS")))
ph <- ph[order(ph$cat != "Spliceosome sig"), ]
ph$sz <- ifelse(ph$cat == "Spliceosome sig", 1.6, 0.4)
p4b <- ggplot(ph, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = cat, size = sz), alpha = 0.55) +
  scale_size_identity() +
  scale_color_manual(values = c("Spliceosome sig" = PURPLE, "Other sig" = "grey80",
                                "Spliceosome NS" = TEAL, "NS" = "grey92")) +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = "logFC (LHS vs LLS)", y = "-log10(p)",
       title = "Phospho: 198/319 (corrected 124/319)")

# (c) kinases
webs2 <- webs
kk <- intersect(colnames(pm), webs2$sample)
kinases <- c("SRPK1","SRPK2","CLK1","CLK2","CLK3","CLK4","CDK7","CDK9","CDK12","CDK13")
kinases <- intersect(kinases, rownames(pm))
kdf <- data.frame(kinase = kinases, rho = NA, p = NA)
for (k in kinases) {
  ct <- suppressWarnings(cor.test(pm[k, kk], webs2$WEBS[match(kk, webs2$sample)], method = "spearman"))
  kdf$rho[kdf$kinase == k] <- ct$estimate; kdf$p[kdf$kinase == k] <- ct$p.value
}
kdf$sig <- ifelse(kdf$p < 0.05, "p<0.05", "n.s.")
p4c <- ggplot(kdf, aes(reorder(kinase, rho), rho, fill = sig)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("p<0.05" = PURPLE, "n.s." = "grey70")) +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = NULL, y = "Spearman rho vs WEBS", title = "Splicing kinases")

# (d) virtual KO (absolute sensitivity ratio)
ko <- read.csv("results/virtual_ko_abs_summary.csv")
ko$type <- ifelse(ko$ko %in% c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP"), "Writer", "Eraser")
p4d <- ggplot(ko, aes(reorder(ko, ratio), ratio, fill = type)) +
  geom_bar(stat = "identity") + coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("Writer" = PURPLE, "Eraser" = TEAL)) +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = NULL, y = "Relative |shift| (spliceosome / other)",
       title = "Virtual KO: spliceosome sensitivity")

fig4 <- (p4a | p4b) / (p4c | p4d) + plot_annotation(tag_levels = "a")
ggsave("figures_GND/Figure_4.pdf", fig4, width = 12, height = 9)
ggsave("figures_GND/Figure_4.png", fig4, width = 12, height = 9, dpi = 200)

# ============ Figure 5: splicing output ============
res <- read.csv("results/splicing_webs_corr.csv")
res$sig_flag <- ifelse(res$q < 0.05 & abs(res$rho) > 0.2, "Sig", "NS")
nsig <- sum(res$sig_flag == "Sig")
p5a <- ggplot(res, aes(rho, -log10(p), color = sig_flag)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Sig" = PURPLE, "NS" = "grey80")) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey40") +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = "Spearman rho (PSI vs WEBS)", y = "-log10(p)",
       title = paste0(nsig, " WEBS-associated AS events (GBM)"))

sig <- res[res$sig_flag == "Sig", ]
tt5 <- as.data.frame(table(sig$splice_type)); colnames(tt5) <- c("type", "count")
p5b <- ggplot(tt5, aes(reorder(type, -count), count, fill = type)) +
  geom_bar(stat = "identity") + scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = "Event type", y = "Count", title = "Event types")

ego <- read.csv("results/splicing_enrichment.csv")
ego <- ego[order(ego$p.adjust), ][1:6, ]
ego$Description <- factor(ego$Description, levels = rev(ego$Description))
p5c <- ggplot(ego, aes(-log10(p.adjust), Description, size = Count, color = p.adjust)) +
  geom_point() + scale_color_gradient(low = PURPLE, high = "grey60") +
  theme_bw(base_size = 9) +
  theme(legend.key.size = unit(0.3, "cm"), axis.text.y = element_text(size = 7)) +
  labs(x = "-log10(adj. p)", y = NULL, size = "Genes", color = "adj. p",
       title = "GO: splicing cascade")

rep <- read.csv("results/splicing_lgg_replication.csv")
ct <- cor.test(rep$rho_gbm, rep$rho_lgg, method = "spearman")
p5d <- ggplot(rep, aes(rho_gbm, rho_lgg)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(alpha = 0.25, size = 0.7, color = PURPLE) +
  theme_bw(base_size = 10) +
  labs(x = "rho in TCGA-GBM", y = "rho in TCGA-LGG",
       title = paste0("LGG replication: 83% concordant, rho = ", round(ct$estimate, 2)))

fig5 <- (p5a | p5b) / (p5c | p5d) + plot_annotation(tag_levels = "a")
ggsave("figures_GND/Figure_5.pdf", fig5, width = 11, height = 9)
ggsave("figures_GND/Figure_5.png", fig5, width = 11, height = 9, dpi = 200)

# ============ Figure 6: clinical validation (4 KM panels) ============
PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)
km_plot <- function(risk, time, event, title) {
  d <- data.frame(risk = risk, time = as.numeric(time), event = as.numeric(event))
  d <- d[!is.na(d$time) & !is.na(d$event) & d$time > 0, ]
  d$grp <- factor(ifelse(d$risk > median(d$risk), "High", "Low"), levels = c("Low", "High"))
  fit <- survfit(Surv(time, event) ~ grp, data = d)
  ggsurvplot(fit, data = d, pval = TRUE, legend.labs = c("Low", "High"),
             palette = c(BLUE, RED), legend = "top",
             title = title, xlab = "Months", ylab = "OS",
             font.main = 11, font.x = 10, font.y = 10)$plot
}
score_cohort <- function(m) {
  g <- intersect(names(PAPER_COEFS), rownames(m))
  as.vector(t(m[g, , drop = FALSE]) %*% PAPER_COEFS[g])
}
e0 <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_expr_log2rsem.rds")
c0 <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")
k1 <- km_plot(score_cohort(as.matrix(e0)), c0$OS.time, c0$OS, "TCGA-GBM (training)")
e693 <- readRDS("data/raw/CGGA/cgga693_expr.rds"); c693 <- readRDS("data/raw/CGGA/cgga693_clinical.rds")
k2 <- km_plot(score_cohort(as.matrix(e693)), c693$OS.time, c693$OS, "CGGA-693")
e325 <- readRDS("data/raw/CGGA/cgga325_expr.rds"); c325 <- readRDS("data/raw/CGGA/cgga325_clinical.rds")
k3 <- km_plot(score_cohort(as.matrix(e325)), c325$OS.time, c325$OS, "CGGA-325")
g1 <- read.csv("results/glass_primary_scores.csv")
k4 <- km_plot(g1$risk, g1$OS.time, g1$OS, "GLASS primary")

fig6 <- (k1 | k2) / (k3 | k4) + plot_annotation(tag_levels = "a")
ggsave("figures_GND/Figure_6.pdf", fig6, width = 11, height = 9)
ggsave("figures_GND/Figure_6.png", fig6, width = 11, height = 9, dpi = 200)

# ============ Figure 7: GLASS recurrence + pharmaco ============
lp <- read.csv("results/glass_recurrence_interval.csv")
ctl <- cor.test(lp$WEBS_TP, lp$surgical_interval_mo, method = "spearman")
p7a <- ggplot(lp, aes(WEBS_TP, surgical_interval_mo)) +
  geom_point(alpha = 0.5, color = PURPLE, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  theme_bw(base_size = 10) +
  labs(x = "Primary tumor WEBS", y = "Months to recurrence",
       title = paste0("GLASS recurrence interval (rho=", round(ctl$estimate, 2), ")"))

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

# pharmaco: GDSC2 HDAC sensitivity + docking summary
gd <- read.csv("results/gdsc2_hdac_sensitivity.csv")
cat("gdsc2 cols:", paste(colnames(gd), collapse = ", "), "\n")
# fallback generic: use gdsc2_drug_webs_corr top hits
gd2 <- read.csv("results/gdsc2_drug_webs_corr.csv")
cat("gdsc2 corr cols:", paste(colnames(gd2), collapse = ", "), "\n")

fig7_top <- (p7a | p7b)
ggsave("figures_GND/Figure_7_partial.pdf", fig7_top, width = 10, height = 4.5)
cat("F3-F7 partial done.\n")
