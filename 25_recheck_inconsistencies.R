# =============================================================================
# Recheck: the 4 inconsistencies vs paper
# =============================================================================
# 1. TCGA-GBM survival stratification (paper: significant; ours: n.s.)
# 2. MLLT1/MED1 dependency (paper: significant; ours: MED1 yes, MLLT1 n.s.)
# 3. Hypoxia decoupling (paper: decoupled; ours: Hallmark rho = -0.45)
# 4. MIF direction (paper: MIF macrophage->malignant axis; ours: up in LLS)
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(Seurat); library(survival); library(survminer)
  library(dplyr); library(ggplot2)
})
set.seed(42)
sink("results/recheck_summary.txt", split = TRUE)

cat("################## 1. TCGA-GBM survival ##################\n\n")

## --- 1a. LinkedOmics legacy cohort (n=144), quartile KM ---
expr <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_expr_log2rsem.rds")
clin <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")

WRITERS <- c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP")
ERASERS <- c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")
webs <- colMeans(expr[intersect(WRITERS, rownames(expr)), ]) -
        colMeans(expr[intersect(ERASERS, rownames(expr)), ])
clin$WEBS <- webs[clin$sample]

q <- quantile(clin$WEBS, c(0.25, 0.75))
clin$grp_q <- ifelse(clin$WEBS <= q[1], "LLS", ifelse(clin$WEBS >= q[2], "LHS", "LDTS"))
clin$grp_med <- ifelse(clin$WEBS > median(clin$WEBS), "High", "Low")

cat("-- LinkedOmics TCGA-GBM (n =", nrow(clin), ") --\n")
cat("Quartile KM (LLS vs LDTS vs LHS):\n")
print(surv_pvalue(survfit(Surv(OS.time, OS) ~ grp_q, data = clin), data = clin))
cat("Median-split KM:\n")
print(surv_pvalue(survfit(Surv(OS.time, OS) ~ grp_med, data = clin), data = clin))
cat("Cox continuous WEBS: ")
print(summary(coxph(Surv(OS.time, OS) ~ WEBS, data = clin))$coefficients)
cat("WEBS range:", range(clin$WEBS), " sd:", sd(clin$WEBS), "\n")

## --- 1b. GDC cohort (STAR counts TPM, n=371 primary) ---
cache <- "data/raw/TCGA-GBM/TCGA-GBM_prepared.rds"
if (file.exists(cache)) {
  dat <- readRDS(cache)
  clin2 <- dat$clinical
  w2 <- intersect(WRITERS, rownames(dat$tpm)); e2 <- intersect(ERASERS, rownames(dat$tpm))
  webs2 <- colMeans(log2(dat$tpm[w2, ] + 1)) - colMeans(log2(dat$tpm[e2, ] + 1))
  clin2$WEBS <- webs2[match(clin2$sample, colnames(dat$tpm))]
  clin2 <- clin2[!is.na(clin2$WEBS), ]
  q2 <- quantile(clin2$WEBS, c(0.25, 0.75), na.rm = TRUE)
  clin2$grp_q <- ifelse(clin2$WEBS <= q2[1], "LLS", ifelse(clin2$WEBS >= q2[2], "LHS", "LDTS"))
  cat("\n-- GDC TCGA-GBM primary (n =", nrow(clin2), ") --\n")
  print(surv_pvalue(survfit(Surv(OS.time, OS) ~ grp_q, data = clin2), data = clin2))
  cat("Cox continuous WEBS:\n")
  print(summary(coxph(Surv(OS.time, OS) ~ WEBS, data = clin2))$coefficients)
}

cat("\n################## 2. MLLT1 / MED1 dependency ##################\n\n")
crispr <- read.csv("data/raw/DepMap/CRISPRGeneEffect.csv", row.names = 1, check.names = FALSE)
colnames(crispr) <- sub(" .*$", "", colnames(crispr))
webs_cl <- read.csv("results/depmap_cellline_webs.csv")
models <- read.csv("data/raw/DepMap/Model.csv", stringsAsFactors = FALSE)

for (set_name in c("all_glioma", "idhwt_gbm")) {
  if (set_name == "all_glioma") {
    lines <- models$ModelID[grepl("glioblastoma|glioma|astrocytoma|gliosarcoma",
                                  models$OncotreePrimaryDisease, ignore.case = TRUE) |
                            grepl("glioblastoma|gliosarcoma", models$OncotreeSubtype,
                                  ignore.case = TRUE)]
  } else {
    lines <- models$ModelID[grepl("glioblastoma|gliosarcoma", models$OncotreeSubtype,
                                  ignore.case = TRUE)]
  }
  use <- intersect(webs_cl$sample, intersect(rownames(crispr), lines))
  wv <- webs_cl$WEBS[match(use, webs_cl$sample)]
  cat("--", set_name, ": n =", length(use), "lines --\n")
  for (g in c("MED1", "MLLT1")) {
    dep <- crispr[use, g]
    ok <- !is.na(dep)
    ct_s <- suppressWarnings(cor.test(wv[ok], dep[ok], method = "spearman"))
    ct_p <- suppressWarnings(cor.test(wv[ok], dep[ok], method = "pearson"))
    cat(sprintf("  %s: spearman r=%.3f p=%.4g | pearson r=%.3f p=%.4g | frac essential (Chronos<-1)=%.2f\n",
                g, ct_s$estimate, ct_s$p.value, ct_p$estimate, ct_p$p.value,
                mean(dep < -1, na.rm = TRUE)))
  }
}

cat("\n################## 3. Hypoxia decoupling (paper Methods exact) ##################\n\n")
obj <- readRDS("data/processed/GSE138794_final.rds")
# Paper: logistic regression Margin-like ~ WEBS + Hallmark_Hypoxia
#        partial correlation: residuals of WEBS ~ Hallmark vs region indicator
md <- obj@meta.data
md$is_margin <- as.integer(md$Region_Type == "Margin-like")
md$WEBS <- as.numeric(md$WEBS)
lr <- glm(is_margin ~ WEBS + Hallmark_Hypoxia, data = md, family = binomial)
cat("Logistic regression Margin ~ WEBS + Hallmark_Hypoxia:\n")
print(summary(lr)$coefficients)

resid_webs <- residuals(lm(WEBS ~ Hallmark_Hypoxia, data = md))
pc <- cor.test(resid_webs, md$is_margin, method = "pearson")
cat("\nPartial correlation (WEBS | hypoxia) vs Margin indicator: r =",
    round(pc$estimate, 3), " p =", signif(pc$p.value, 3), "\n")

sp_all <- cor.test(md$WEBS, md$Hallmark_Hypoxia, method = "spearman")
cat("Raw Spearman WEBS vs Hallmark Hypoxia: rho =", round(sp_all$estimate, 3),
    " p =", signif(sp_all$p.value, 3), "\n")
# within malignant only
mal <- md[obj$is_malignant == TRUE & !is.na(obj$is_malignant), ]
if (nrow(mal) > 100) {
  sp_m <- cor.test(mal$WEBS, mal$Hallmark_Hypoxia, method = "spearman")
  cat("Malignant-only Spearman WEBS vs Hallmark: rho =", round(sp_m$estimate, 3),
      " p =", signif(sp_m$p.value, 3), "\n")
}

cat("\n################## 4. MIF axis (paper Methods exact) ##################\n\n")
# interaction score = mean(MIF in macrophage/microglia) x mean(CD74+CXCR4 in malignant)
# permutation test: 10,000 label shuffles; compare across WEBS groups
data_mat <- GetAssayData(obj, layer = "data")
is_mac <- (!is.na(obj$original_type) & obj$original_type == "Myeloid") |
          (is.na(obj$original_type) & obj$cell_type == "Macrophage_Microglia")
is_mal <- (!is.na(obj$original_type) & obj$original_type == "Neoplastic_cell") |
          (is.na(obj$original_type) & obj$cell_type == "Malignant")
if ("MIF" %in% rownames(data_mat)) {
  mif_mac <- mean(data_mat["MIF", is_mac])
  cd74_mal <- mean(data_mat["CD74", is_mal])
  cxcr4_mal <- if ("CXCR4" %in% rownames(data_mat)) mean(data_mat["CXCR4", is_mal]) else NA
  cat("MIF (myeloid mean) =", round(mif_mac, 3),
      " | CD74 (malignant mean) =", round(cd74_mal, 3),
      " | CXCR4 (malignant mean) =", round(cxcr4_mal, 3), "\n")

  # per WEBS group: interaction score
  for (grp in c("LLS", "LDTS", "LHS")) {
    cells <- WhichCells(obj, expression = WEBS_group == grp)
    mac_g <- intersect(cells, colnames(obj)[is_mac])
    mal_g <- intersect(cells, colnames(obj)[is_mal])
    if (length(mac_g) > 50 && length(mal_g) > 50) {
      score <- mean(data_mat["MIF", mac_g]) *
               mean(colMeans(data_mat[c("CD74", "CXCR4"), mal_g, drop = FALSE]))
      cat(grp, ": MIF->CD74/CXCR4 interaction score =", round(score, 4),
          " (n_mac =", length(mac_g), ", n_mal =", length(mal_g), ")\n")
    }
  }

  # permutation: shuffle cell-type labels, null distribution of the score
  n_perm <- 1000
  null_scores <- replicate(n_perm, {
    mac_s <- sample(colnames(data_mat), sum(is_mac))
    mal_s <- sample(colnames(data_mat), sum(is_mal))
    mean(data_mat["MIF", mac_s]) * mean(colMeans(data_mat[c("CD74","CXCR4"), mal_s, drop = FALSE]))
  })
  obs <- mif_mac * mean(colMeans(data_mat[c("CD74","CXCR4"), is_mal, drop = FALSE]))
  p_perm <- (sum(null_scores >= obs) + 1) / (n_perm + 1)
  cat("\nObserved MIF->CD74/CXCR4 score:", round(obs, 4),
      " | permutation p (null >= obs):", round(p_perm, 4), "\n")
} else {
  cat("MIF gene not found in expression matrix\n")
}

sink()
message("Recheck complete -> results/recheck_summary.txt")
