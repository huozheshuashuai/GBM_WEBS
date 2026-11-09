# =============================================================================
# 79. BJC layout: merge GLASS longitudinal into Figure 6 (a-f);
#     pharmaco panels -> new supplementary figure S35
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork); library(survival); library(survminer) })
PURPLE <- "#7B4EA3"; TEAL <- "#3D9B8C"; RED <- "#d62728"; BLUE <- "#1f77b4"

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
  ggsurvplot(fit, data = d, pval = TRUE, linetype = "strata",
             legend.labs = c(paste0("Low (n=", nlo, ")"), paste0("High (n=", nhi, ")")),
             palette = c(BLUE, RED), legend = "top",
             title = title, subtitle = paste0("HR = ", hr, " (95% CI ", lo, "-", hi, ")"),
             xlab = "Months", ylab = "OS",
             font.main = 10, font.subtitle = 8, font.x = 9, font.y = 9,
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

# GLASS recurrence interval (first-pair)
lp <- read.csv("results/glass_recurrence_interval.csv")
fp <- lp[!duplicated(lp$case_barcode), ]
ctf <- cor.test(fp$WEBS_TP, fp$surgical_interval_mo, method = "spearman")
p6e <- ggplot(fp, aes(WEBS_TP, surgical_interval_mo)) +
  geom_point(alpha = 0.6, color = PURPLE, size = 1.2) +
  geom_smooth(method = "loess", se = TRUE, color = "grey30", formula = y ~ x) +
  theme_bw(base_size = 9) +
  labs(x = "Primary tumor WEBS", y = "Months to first recurrence",
       title = "GLASS: recurrence interval (exploratory)",
       subtitle = paste0("one record/patient (n = ", nrow(fp), "); rho = ", round(ctf$estimate, 2)))

# GLASS longitudinal pairs
lo <- read.csv("results/glass_longitudinal_pairs.csv")
plot_df <- data.frame(case = rep(lo$case_barcode, 2),
                      stage = factor(rep(c("Primary", "Recurrence"), each = nrow(lo)),
                                     levels = c("Primary", "Recurrence")),
                      WEBS = c(lo$WEBS_TP, lo$WEBS_R1))
p6f <- ggplot(plot_df, aes(stage, WEBS, group = case)) +
  geom_line(alpha = 0.25, color = "grey50") +
  geom_point(aes(color = stage), size = 1) +
  scale_color_manual(values = c("Primary" = TEAL, "Recurrence" = PURPLE), name = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom") +
  labs(x = NULL, y = "WEBS",
       title = paste0("GLASS pairs (n = ", nrow(lo), " pairs / ",
                      length(unique(lo$case_barcode)), " patients)"))

fig6 <- ((k1$plot / k1$table) | (k2$plot / k2$table) | (k3$plot / k3$table)) /
        ((k4$plot / k4$table) | p6e | p6f) +
        plot_annotation(tag_levels = "a")
ggsave("figures_BJC/Figure_6.pdf", fig6, width = 15, height = 11)
ggsave("figures_BJC/Figure_6.png", fig6, width = 15, height = 11, dpi = 200)

# ---- S35: pharmaco (old Fig 7c + 7d) ----
gd <- read.csv("results/gdsc2_drug_webs_corr.csv")
gd <- gd[!is.na(gd$q) & gd$q < 0.1, ]
top_drugs <- gd[order(gd$spearman_r), ][1:10, ]
top_drugs$DRUG_NAME <- factor(top_drugs$DRUG_NAME, levels = rev(top_drugs$DRUG_NAME))
p35a <- ggplot(top_drugs, aes(DRUG_NAME, spearman_r)) +
  geom_bar(stat = "identity", fill = PURPLE) + coord_flip() +
  theme_bw(base_size = 10) +
  labs(x = NULL, y = "Spearman rho (WEBS vs ln IC50)",
       title = "Top WEBS-sensitizing single agents",
       subtitle = "GDSC2, q < 0.1; pan-cancer cell-line pharmacogenomic associations",
       caption = "Single-agent associations; combination predictions are separate")

dd <- read.csv("results/depmap_webs_dependency.csv")
dd$gene <- dd[[1]]
key <- dd[dd$gene %in% c("MED1","SF3B1","SRSF1","HDAC2","SIRT2"), ]
key$sig <- ifelse(key$q < 0.05, "FDR q < 0.05",
                  ifelse(key$p < 0.05, "nominal p < 0.05", "n.s."))
key$gene <- factor(key$gene, levels = key$gene[order(key$spearman_r)])
p35b <- ggplot(key, aes(gene, spearman_r, fill = sig)) +
  geom_bar(stat = "identity") + coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("FDR q < 0.05" = PURPLE, "nominal p < 0.05" = "#D4A017",
                               "n.s." = "grey65"), name = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom") +
  labs(x = NULL, y = "Correlation with dependency (Chronos)",
       title = "DepMap dependency (glioma lines)")

fig_s35 <- (p35a | p35b) + plot_annotation(tag_levels = "a")
ggsave("figures_BJC/FigS35_pharmaco.pdf", fig_s35, width = 11, height = 5)
cat("Fig6 merged + S35 built.\n")
