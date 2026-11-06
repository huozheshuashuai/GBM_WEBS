# =============================================================================
# 39. My four original analyses
#  (1) Jonckheere trend test: acetylation sites monotone LLS<LDTS<LHS
#  (2) Proteome-wide discordance tail: spliceosome vs all proteins
#  (3) Protein half-life cross-validation (Mathieson 2018, if downloadable)
#  (4) Random-effects meta-analysis of the six-gene signature across cohorts
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(dplyr); library(ggplot2); library(clinfun); library(meta)
})
set.seed(42)
dir.create("results", showWarnings = FALSE)
dir.create("figures_final", showWarnings = FALSE)

sink("results/my_analyses_summary.txt", split = TRUE)

webs <- read.csv("results/cptac_webs.csv")

# ============ (1) Jonckheere trend test: LLS < LDTS < LHS ============
acetyl <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_acetylome_mssm_site.cct",
                               row.names = 1, check.names = FALSE))
grp <- webs$WEBS_group[match(colnames(acetyl), webs$sample)]
ok <- !is.na(grp)
acetyl <- acetyl[, ok]; grp <- grp[ok]
grp_ord <- factor(grp, levels = c("LLS", "LDTS", "LHS"), ordered = TRUE)

de <- read.csv("results/acetyl_sites_limma.csv")
sig_splice <- de %>% filter(is_spliceosome, !is.na(adj.P.Val), adj.P.Val < 0.05) %>% pull(site)
sig_splice <- intersect(sig_splice, rownames(acetyl))

jt_res <- lapply(sig_splice, function(s) {
  x <- acetyl[s, ]; keep <- !is.na(x)
  if (sum(keep) < 15) return(NULL)
  jt <- jonckheere.test(x[keep], grp_ord[keep], alternative = "increasing")
  data.frame(site = s, p = jt$p.value)
})
jt_df <- do.call(rbind, jt_res)
jt_df$q <- p.adjust(jt_df$p, method = "BH")
write.csv(jt_df, "results/jonckheere_trend.csv", row.names = FALSE)
cat("=== (1) Jonckheere trend test (LLS<LDTS<LHS) ===\n")
cat("Significant spliceosome sites tested:", nrow(jt_df), "\n")
cat("Monotone increasing (FDR<0.05):", sum(jt_df$q < 0.05), "\n")
cat("Median p:", signif(median(jt_df$p), 3), "\n\n")

# ============ (2) Proteome-wide discordance tail ============
proteome <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct",
                                 row.names = 1, check.names = FALSE))
rnaseq <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct",
                               row.names = 1, check.names = FALSE))
source("tools/WEBS/webs.R")
source("08_cptac_proteomics.R")
common <- intersect(colnames(proteome), intersect(colnames(rnaseq), webs$sample))
prot_w <- webs$WEBS[match(common, webs$sample)]
rna_webs <- webs_compute_bulk(rnaseq)$WEBS[match(common, webs$sample)]

all_genes <- intersect(rownames(proteome), rownames(rnaseq))
all_genes <- all_genes[all_genes != "" & !is.na(all_genes)]
splice_genes <- intersect(get_spliceosome_genes(), all_genes)

disc <- sapply(all_genes, function(g) {
  okp <- !is.na(proteome[g, common]); okr <- !is.na(rnaseq[g, common])
  if (sum(okp) < 20 || sum(okr) < 20) return(NA)
  rp <- cor(prot_w[okp], proteome[g, common][okp], method = "spearman")
  rm_ <- cor(rna_webs[okr], rnaseq[g, common][okr], method = "spearman")
  rp - rm_
})
names(disc) <- all_genes
disc <- disc[!is.na(disc)]
is_spl <- names(disc) %in% splice_genes
emp_p <- (sum(disc[!is_spl] >= mean(disc[is_spl])) + 1) / (sum(!is_spl) + 1)
cat("=== (2) Proteome-wide discordance tail ===\n")
cat("Genes tested:", length(disc), "| spliceosome:", sum(is_spl), "\n")
cat("Spliceosome mean discordance:", round(mean(disc[is_spl]), 3),
    " | all others:", round(mean(disc[!is_spl]), 3), "\n")
cat("Wilcoxon p:", signif(wilcox.test(disc[is_spl], disc[!is_spl])$p.value, 3), "\n")
cat("Empirical p (spliceosome mean >= others):", signif(emp_p, 3), "\n\n")

df_plot <- data.frame(gene = names(disc), discordance = disc, spliceosome = is_spl)
p <- ggplot(df_plot, aes(x = spliceosome, y = discordance)) +
  geom_violin(aes(fill = spliceosome), scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "#C0392B")) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(x = NULL, y = "Discordance (protein rho - mRNA rho vs WEBS)",
       title = "Spliceosome proteins sit at the discordance extreme",
       subtitle = sprintf("Wilcoxon p = %.2g (n = %d spliceosome vs %d others)",
                          wilcox.test(disc[is_spl], disc[!is_spl])$p.value,
                          sum(is_spl), sum(!is_spl))) +
  scale_x_discrete(labels = c("FALSE" = "All other proteins", "TRUE" = "Spliceosome"))
ggsave("figures_final/fig_discordance_tail.pdf", p, width = 6.5, height = 5.5)

# ============ (3) Protein half-life cross-validation ============
cat("=== (3) Protein half-life (Mathieson 2018) ===\n")
hl_url <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302754-mmc2.xlsx"
hl_file <- "data/raw/mathieson_halflife.xlsx"
if (!file.exists(hl_file)) {
  tryCatch(download.file(hl_url, hl_file, mode = "wb", quiet = TRUE),
           error = function(e) message("download failed"))
}
if (file.exists(hl_file) && file.size(hl_file) > 10000) {
  library(readxl)
  hl <- tryCatch(read_excel(hl_file), error = function(e) NULL)
  if (!is.null(hl)) {
    hlcols <- colnames(hl)
    genecol <- hlcols[grepl("gene", hlcols, ignore.case = TRUE)][1]
    hlcol <- hlcols[grepl("half", hlcols, ignore.case = TRUE)][1]
    hl <- hl[!is.na(hl[[genecol]]), ]
    hl$is_splice <- toupper(hl[[genecol]]) %in% splice_genes
    wt <- wilcox.test(as.numeric(hl[[hlcol]][hl$is_splice]),
                      as.numeric(hl[[hlcol]][!hl$is_splice]))
    cat("Spliceosome median half-life:",
        round(median(as.numeric(hl[[hlcol]][hl$is_splice]), na.rm = TRUE), 1), "\n")
    cat("Others median half-life:",
        round(median(as.numeric(hl[[hlcol]][!hl$is_splice]), na.rm = TRUE), 1), "\n")
    cat("Wilcoxon p:", signif(wt$p.value, 3), "\n")
    write.csv(data.frame(gene = hl[[genecol]], half_life = hl[[hlcol]],
                         spliceosome = hl$is_splice),
              "results/halflife_crossval.csv", row.names = FALSE)
  } else cat("file unreadable\n")
} else cat("half-life data unavailable (download failed) — skipped\n")

# ============ (4) Random-effects meta-analysis ============
cat("\n=== (4) Meta-analysis of six-gene signature ===\n")
PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)
library(survival)
cox_cohort <- function(expr, clin, cohort) {
  g <- intersect(names(PAPER_COEFS), rownames(expr))
  rs <- as.vector(t(expr[g, , drop = FALSE]) %*% PAPER_COEFS[g])
  clin$risk <- rs[match(clin$sample, colnames(expr))]
  clin <- clin[!is.na(clin$risk), ]
  # scale risk for comparable HR
  clin$risk <- as.numeric(scale(clin$risk))
  m <- coxph(Surv(OS.time, OS) ~ risk, data = clin)
  s <- summary(m)
  data.frame(cohort = cohort, TE = coef(m)["risk"],
             seTE = s$coefficients["risk", 3])
}
res <- list()
e <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_expr_log2rsem.rds")
c <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")
res[[1]] <- cox_cohort(e, c, "TCGA-GBM")
for (cg in c("cgga693", "cgga325")) {
  res[[length(res) + 1]] <- cox_cohort(readRDS(paste0("data/raw/CGGA/", cg, "_expr.rds")),
                                       readRDS(paste0("data/raw/CGGA/", cg, "_clinical.rds")), cg)
}
res <- do.call(rbind, res)
write.csv(res, "results/meta_analysis_input.csv", row.names = FALSE)

m <- metagen(TE = res$TE, seTE = res$seTE, studlab = res$cohort,
             sm = "HR", method.tau = "REML")
cat("\nRandom-effects pooled HR (per SD of risk score):\n")
print(summary(m))
pdf("figures_final/fig_meta_forest.pdf", width = 7, height = 4)
forest(m, xlab = "Hazard Ratio (per SD risk score)",
       leftlabs = c("Cohort"), rightlabs = c("HR [95% CI]"))
dev.off()

sink()
cat("\nDone. See results/my_analyses_summary.txt\n")
