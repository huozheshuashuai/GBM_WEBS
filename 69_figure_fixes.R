# =============================================================================
# 69. GND figure fixes round: F4 (legends + GENIE3 naming), F5 (layout + titles),
#     F6 (risk tables + HR/CI + n), F7 (first-pair GLASS + legends)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork); library(survival); library(survminer); library(tidyr) })
PURPLE <- "#7B4EA3"; TEAL <- "#3D9B8C"; RED <- "#d62728"; BLUE <- "#1f77b4"
dir.create("figures_GND", showWarnings = FALSE)
source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
webs <- read.csv("results/cptac_webs.csv")
prot <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct", check.names = FALSE)
pm <- as.matrix(sapply(prot[, -1], as.numeric)); rownames(pm) <- prot[[1]]

# ============ Figure 4 (fixed) ============
tt <- read.csv("results/acetyl_sites_limma.csv")
tt$gene <- sub("-.*$", "", tt$site)
top <- tt %>% filter(gene %in% splice, adj.P.Val < 0.05) %>% arrange(desc(abs(t))) %>% head(30)
ac <- read.delim("E:/KIMI/data/HS_CPTAC_GBM_acetylome_mssm_site.cct", check.names = FALSE)
sub <- ac[match(top$site, ac[[1]]), ]
am <- as.matrix(sapply(sub[, -1], as.numeric)); rownames(am) <- top$site
common <- intersect(colnames(am), webs$sample)
am <- am[, common]; w2 <- webs[match(common, webs$sample), ]
am <- am[, order(w2$WEBS)]
grp_ord <- w2$WEBS_group[order(w2$WEBS)]
long <- pivot_longer(as.data.frame(am) %>% mutate(site = rownames(am)), -site,
                     names_to = "sample", values_to = "val")
long$site <- factor(long$site, levels = rev(top$site))
long$sample <- factor(long$sample, levels = colnames(am))
# annotation strip data
ann <- data.frame(sample = colnames(am), group = grp_ord)
ann$sample <- factor(ann$sample, levels = colnames(am))
p_ann <- ggplot(ann, aes(sample, 1, fill = group)) +
  geom_tile() +
  scale_fill_manual(values = c("LLS" = BLUE, "LDTS" = "#ff7f0e", "LHS" = RED), name = "WEBS group") +
  theme_void(base_size = 9) +
  theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"))
p4a_main <- ggplot(long, aes(sample, site, fill = val)) +
  geom_tile() +
  scale_fill_gradient2(low = TEAL, mid = "white", high = PURPLE, midpoint = 0,
                       na.value = "grey90", name = "Acetyl.") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6.5)) +
  labs(x = NULL, y = NULL, title = "Spliceosome acetylation (OR = 7.36; corrected 2.29)")
p4a <- (p4a_main / p_ann + plot_layout(heights = c(20, 1))) + plot_annotation(tag_levels = list(c("", "")))

ph <- read.csv("results/phospho_limma_LHS_vs_LLS.csv")
ph$cat <- ifelse(ph$adj.P.Val < 0.05 & ph$is_spliceosome & ph$logFC > 0, "Spliceosome sig (up)",
          ifelse(ph$adj.P.Val < 0.05 & ph$is_spliceosome, "Spliceosome sig (down)",
          ifelse(ph$adj.P.Val < 0.05, "Other sig",
          ifelse(ph$is_spliceosome, "Spliceosome NS", "NS"))))
ph <- ph[order(!grepl("Spliceosome sig", ph$cat)), ]
ph$sz <- ifelse(grepl("Spliceosome sig", ph$cat), 1.6, 0.4)
p4b <- ggplot(ph, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = cat, size = sz), alpha = 0.55) +
  scale_size_identity() +
  scale_color_manual(values = c("Spliceosome sig (up)" = PURPLE, "Spliceosome sig (down)" = "#1b9e77",
                                "Other sig" = "grey80", "Spliceosome NS" = TEAL, "NS" = "grey92"),
                     name = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7.5)) +
  guides(color = guide_legend(override.aes = list(size = 2), nrow = 2)) +
  labs(x = "logFC (high- vs low-WEBS)", y = "-log10(p)",
       title = "Phosphosites: 198/319 sig (corrected 124/319)")

kk <- intersect(colnames(pm), webs$sample)
kinases <- intersect(c("SRPK1","SRPK2","CLK1","CLK2","CLK3","CLK4","CDK7","CDK9","CDK12","CDK13"), rownames(pm))
kdf <- data.frame(kinase = kinases, rho = NA, p = NA)
for (k in kinases) {
  ct <- suppressWarnings(cor.test(pm[k, kk], webs$WEBS[match(kk, webs$sample)], method = "spearman"))
  kdf$rho[kdf$kinase == k] <- ct$estimate; kdf$p[kdf$kinase == k] <- ct$p.value
}
kdf$sig <- ifelse(kdf$p < 0.05, "p < 0.05", "n.s.")
p4c <- ggplot(kdf, aes(reorder(kinase, rho), rho, fill = sig)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("p < 0.05" = PURPLE, "n.s." = "grey70"), name = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom") +
  labs(x = NULL, y = "Spearman rho vs WEBS", title = "Splicing-kinase abundance vs WEBS")

ko <- read.csv("results/virtual_ko_abs_summary.csv")
ko$type <- ifelse(ko$ko %in% c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP"), "Writers", "Erasers")
p4d <- ggplot(ko, aes(reorder(ko, ratio), ratio, fill = type)) +
  geom_bar(stat = "identity") + coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("Writers" = PURPLE, "Erasers" = TEAL), name = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom") +
  labs(x = NULL, y = "Relative |shift| (spliceosome / other)",
       title = "GENIE3 connectivity sensitivity")

fig4 <- ((p4a | p4b) / (p4c | p4d)) + plot_annotation(tag_levels = list(c("a","","b","c","d")))
ggsave("figures_GND/Figure_4.pdf", fig4, width = 12, height = 10)
ggsave("figures_GND/Figure_4.png", fig4, width = 12, height = 10, dpi = 200)

# ============ Figure 5 (fixed) ============
res <- read.csv("results/splicing_webs_corr.csv")
res$sig_flag <- ifelse(res$q < 0.05 & abs(res$rho) > 0.2, "Significant", "NS")
nsig <- sum(res$sig_flag == "Significant")
xmax <- max(abs(res$rho), na.rm = TRUE) * 1.05
p5a <- ggplot(res, aes(rho, -log10(p), color = sig_flag)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Significant" = PURPLE, "NS" = "grey80"), name = NULL) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey40") +
  scale_x_continuous(limits = c(-xmax, xmax)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", plot.margin = margin(5, 15, 5, 5)) +
  labs(x = "Spearman rho (PSI vs WEBS)", y = "-log10(p)",
       title = paste0(nsig, " WEBS-associated AS events"),
       caption = "Dashed lines: |rho| = 0.2 selection threshold")

sig <- res[res$sig_flag == "Significant", ]
tt5 <- as.data.frame(table(sig$splice_type)); colnames(tt5) <- c("type", "count")
p5b <- ggplot(tt5, aes(reorder(type, -count), count, fill = type)) +
  geom_bar(stat = "identity") + scale_fill_brewer(palette = "Set2", name = "Event type") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7.5)) +
  labs(x = "Event type", y = "Count", title = "Event-type composition")

ego <- read.csv("results/splicing_enrichment.csv")
ego <- ego[order(ego$p.adjust), ][1:6, ]
ego$Description <- factor(ego$Description, levels = rev(ego$Description))
p5c <- ggplot(ego, aes(-log10(p.adjust), Description, size = Count, color = p.adjust)) +
  geom_point() + scale_color_gradient(low = PURPLE, high = "grey60") +
  theme_bw(base_size = 9) +
  theme(legend.key.size = unit(0.3, "cm"), axis.text.y = element_text(size = 7),
        legend.position = "right") +
  labs(x = "-log10(adj. p)", y = NULL, size = "Genes", color = "adj. p",
       title = "GO enrichment of genes harboring\nWEBS-associated splicing events")

rep <- read.csv("results/splicing_lgg_replication.csv")
ctr <- cor.test(rep$rho_gbm, rep$rho_lgg, method = "spearman")
p5d <- ggplot(rep, aes(rho_gbm, rho_lgg)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(alpha = 0.25, size = 0.7, color = PURPLE) +
  theme_bw(base_size = 10) +
  labs(x = "rho in TCGA-GBM", y = "rho in TCGA-LGG",
       title = "Cross-cohort replication (TCGA-LGG)",
       subtitle = paste0("83% sign-concordant; effect-size rho = ", round(ctr$estimate, 2)))

fig5 <- (p5a | p5b) / (p5c | p5d) + plot_annotation(tag_levels = "a")
ggsave("figures_GND/Figure_5.pdf", fig5, width = 11, height = 10)
ggsave("figures_GND/Figure_5.png", fig5, width = 11, height = 10, dpi = 200)

# ============ Figure 6 (fixed: risk tables + n + HR) ============
PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)
score_cohort <- function(m) {
  g <- intersect(names(PAPER_COEFS), rownames(m))
  as.vector(t(m[g, , drop = FALSE]) %*% PAPER_COEFS[g])
}
km_plot2 <- function(risk, time, event, title) {
  d <- data.frame(risk = risk, time = as.numeric(time), event = as.numeric(event))
  d <- d[!is.na(d$time) & !is.na(d$event) & d$time > 0, ]
  d$grp <- factor(ifelse(d$risk > median(d$risk), "High", "Low"), levels = c("Low", "High"))
  fit <- survfit(Surv(time, event) ~ grp, data = d)
  cox <- summary(coxph(Surv(time, event) ~ grp, data = d))
  hr <- round(cox$conf.int[1, 1], 2); lo <- round(cox$conf.int[1, 3], 2); hi <- round(cox$conf.int[1, 4], 2)
  nlo <- sum(d$grp == "Low"); nhi <- sum(d$grp == "High")
  ggsurvplot(fit, data = d, pval = TRUE, pval.method = TRUE,
             legend.labs = c(paste0("Low (n=", nlo, ")"), paste0("High (n=", nhi, ")")),
             palette = c(BLUE, RED), legend = "top",
             title = title, subtitle = paste0("HR = ", hr, " (95% CI ", lo, "-", hi, ")"),
             xlab = "Months", ylab = "Overall survival",
             font.main = 11, font.subtitle = 9, font.x = 10, font.y = 10,
             risk.table = TRUE, risk.table.height = 0.22,
             risk.table.title = "Number at risk", tables.theme = theme_cleantable())
}
e0 <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_expr_log2rsem.rds")
c0 <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")
k1 <- km_plot2(score_cohort(as.matrix(e0)), c0$OS.time, c0$OS, "TCGA-GBM (training)")
e693 <- readRDS("data/raw/CGGA/cgga693_expr.rds"); c693 <- readRDS("data/raw/CGGA/cgga693_clinical.rds")
k2 <- km_plot2(score_cohort(as.matrix(e693)), c693$OS.time, c693$OS, "CGGA-693")
e325 <- readRDS("data/raw/CGGA/cgga325_expr.rds"); c325 <- readRDS("data/raw/CGGA/cgga325_clinical.rds")
k3 <- km_plot2(score_cohort(as.matrix(e325)), c325$OS.time, c325$OS, "CGGA-325")
g1 <- read.csv("results/glass_primary_scores.csv")
k4 <- km_plot2(g1$risk, g1$OS.time, g1$OS, "GLASS primary")

fig6 <- ((k1$plot / k1$table) | (k2$plot / k2$table)) /
        ((k3$plot / k3$table) | (k4$plot / k4$table)) + plot_annotation(tag_levels = list(c("a","","b","","c","","d","")))
ggsave("figures_GND/Figure_6.pdf", fig6, width = 12, height = 11)
ggsave("figures_GND/Figure_6.png", fig6, width = 12, height = 11, dpi = 200)

# ============ Figure 7 (fixed: first-pair + legends) ============
lp <- read.csv("results/glass_recurrence_interval.csv")
fp <- lp[!duplicated(lp$case_barcode), ]
ctf <- cor.test(fp$WEBS_TP, fp$surgical_interval_mo, method = "spearman")
p7a <- ggplot(fp, aes(WEBS_TP, surgical_interval_mo)) +
  geom_point(alpha = 0.6, color = PURPLE, size = 1.4) +
  geom_smooth(method = "loess", se = TRUE, color = "grey30", formula = y ~ x) +
  theme_bw(base_size = 10) +
  labs(x = "Primary tumor WEBS", y = "Months to first recurrence",
       title = "GLASS: recurrence interval",
       subtitle = paste0("one pair per patient (n = ", nrow(fp), "); Spearman rho = ",
                         round(ctf$estimate, 2), ", p = ", signif(ctf$p.value, 2)),
       caption = "Trend line (LOESS) is a visual aid; the test is Spearman correlation")

lo <- read.csv("results/glass_longitudinal_pairs.csv")
n_pat <- length(unique(lo$case_barcode))
plot_df <- data.frame(case = rep(lo$case_barcode, 2),
                      stage = factor(rep(c("Primary", "Recurrence"), each = nrow(lo)),
                                     levels = c("Primary", "Recurrence")),
                      WEBS = c(lo$WEBS_TP, lo$WEBS_R1))
p7b <- ggplot(plot_df, aes(stage, WEBS, group = case)) +
  geom_line(alpha = 0.25, color = "grey50") +
  geom_point(aes(color = stage), size = 1.2) +
  scale_color_manual(values = c("Primary" = TEAL, "Recurrence" = PURPLE), name = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom") +
  labs(x = NULL, y = "WEBS",
       title = "GLASS longitudinal pairs",
       subtitle = paste0(nrow(lo), " pairs from ", n_pat, " patients; paired Wilcoxon p = 0.09"))

gd <- read.csv("results/gdsc2_drug_webs_corr.csv")
gd <- gd[!is.na(gd$q) & gd$q < 0.1, ]
top_drugs <- gd[order(gd$spearman_r), ][1:10, ]
top_drugs$DRUG_NAME <- factor(top_drugs$DRUG_NAME, levels = rev(top_drugs$DRUG_NAME))
p7c <- ggplot(top_drugs, aes(DRUG_NAME, spearman_r)) +
  geom_bar(stat = "identity", fill = PURPLE) + coord_flip() +
  theme_bw(base_size = 10, ) +
  theme(axis.text.y = element_text(size = 9)) +
  labs(x = NULL, y = "Spearman rho (WEBS vs ln IC50)",
       title = "Top WEBS-sensitizing single agents",
       subtitle = "GDSC2, FDR q < 0.1; negative rho = higher sensitivity",
       caption = "Single-agent associations; combination predictions are separate")

dd <- read.csv("results/depmap_webs_dependency.csv")
dd$gene <- dd[[1]]
key <- dd[dd$gene %in% c("MED1","SF3B1","SRSF1","HDAC2","SIRT2"), ]
key$sig <- ifelse(key$q < 0.05, "FDR q < 0.05",
                  ifelse(key$p < 0.05, "nominal p < 0.05", "n.s."))
key$gene <- factor(key$gene, levels = key$gene[order(key$spearman_r)])
p7d <- ggplot(key, aes(gene, spearman_r, fill = sig)) +
  geom_bar(stat = "identity") + coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("FDR q < 0.05" = PURPLE, "nominal p < 0.05" = "#D4A017",
                               "n.s." = "grey65"), name = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom") +
  labs(x = NULL, y = "Correlation with dependency (Chronos)",
       title = "DepMap dependency (glioma lines)",
       subtitle = "Negative = stronger dependency in high-WEBS lines")

fig7 <- (p7a | p7b) / (p7c | p7d) + plot_annotation(tag_levels = "a")
ggsave("figures_GND/Figure_7.pdf", fig7, width = 12, height = 10)
ggsave("figures_GND/Figure_7.png", fig7, width = 12, height = 10, dpi = 200)
cat("F4-F7 fixed.\n")
