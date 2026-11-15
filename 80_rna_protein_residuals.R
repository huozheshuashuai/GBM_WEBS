# =============================================================================
# 80. DIRECT test of post-transcriptional accumulation:
#     per-protein residuals (protein abundance beyond what mRNA predicts)
#     Does the spliceosome show systematically positive residuals in high-WEBS?
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)
source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()

prot <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct", check.names = FALSE)
rna <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct", check.names = FALSE)
pm <- as.matrix(sapply(prot[, -1], as.numeric)); rownames(pm) <- prot[[1]]
rm_ <- as.matrix(sapply(rna[, -1], as.numeric)); rownames(rm_) <- rna[[1]]
webs <- read.csv("results/cptac_webs.csv")

genes <- intersect(rownames(pm), rownames(rm_))
samples <- intersect(colnames(pm), intersect(colnames(rm_), webs$sample))
pm <- pm[genes, samples]; rm_ <- rm_[genes, samples]
w <- webs$WEBS[match(samples, webs$sample)]
grp <- factor(webs$WEBS_group[match(samples, webs$sample)], levels = c("LLS","LDTS","LHS"))

# ---- per-protein: lm(protein ~ mRNA) across patients; residual = protein - predicted ----
# model A per-protein coefficient of WEBS after mRNA adjustment
beta_webs <- res_list <- list()
for (g in genes) {
  fit <- lm(pm[g, ] ~ rm_[g, ])
  res_list[[g]] <- resid(fit)
}
resid_m <- do.call(rbind, res_list)   # proteins x samples, mRNA-adjusted residuals
# WEBS coefficient after mRNA adjustment (model A)
get_beta <- function(g) {
  fit <- tryCatch(lm(pm[g, ] ~ rm_[g, ] + w), error = function(e) NULL)
  if (is.null(fit)) return(c(NA, NA))
  co <- summary(fit)$coefficients
  if (nrow(co) < 3) return(c(NA, NA))
  c(co["w", "Estimate"], co["w", "Pr(>|t|)"])
}
ba <- t(sapply(genes, get_beta))
beta_df <- data.frame(gene = genes, beta_webs = ba[, 1], p = ba[, 2])
beta_df$q <- p.adjust(beta_df$p, "BH")
beta_df$is_spliceosome <- beta_df$gene %in% splice
write.csv(beta_df, "results/rna_protein_residual_perprotein.csv", row.names = FALSE)

# ---- model B: spliceosome-complex residual score per sample ----
sp <- intersect(splice, genes)
splice_resid_score <- colMeans(resid_m[sp, , drop = FALSE], na.rm = TRUE)
grp2 <- grp
df <- data.frame(WEBS = w, group = grp2, resid_score = splice_resid_score)
ct <- cor.test(w, splice_resid_score, method = "spearman")
anova_p <- summary(aov(resid_score ~ group, data = df))[[1]][1, "Pr(>F)"]
kw <- kruskal.test(resid_score ~ group, data = df)

# specificity: abundance- AND variance-matched random sets on residual score
# stricter control: background pool excludes spliceosome genes; sampling is
# without replacement; genes are matched on mean abundance decile x within-decile
# variance tercile (fallback to the whole abundance decile if a bin runs out)
expr_mean <- rowMeans(pm, na.rm = TRUE)
expr_sd <- apply(pm, 1, sd, na.rm = TRUE)
dec <- cut(expr_mean, breaks = quantile(expr_mean, probs = 0:10/10, na.rm = TRUE),
           include.lowest = TRUE, labels = FALSE)
names(dec) <- genes
bin <- integer(length(genes)); names(bin) <- genes
for (d in 1:10) {
  idx <- which(dec == d)
  vt <- cut(expr_sd[idx], breaks = unique(quantile(expr_sd[idx], probs = 0:3/3, na.rm = TRUE)),
            include.lowest = TRUE, labels = FALSE)
  bin[idx] <- (d - 1) * 3 + as.integer(vt)
}
pool <- setdiff(genes, sp)                      # spliceosome excluded from background
sp_bins <- bin[sp]
bin_levels <- sort(unique(sp_bins))
bin_pool <- lapply(bin_levels, function(b) which(bin == b & names(bin) %in% pool))
names(bin_pool) <- bin_levels
dec_pool <- lapply(1:10, function(d) which(dec == d & names(dec) %in% pool))
bin_dec <- ceiling(bin_levels / 3)
n_need <- sapply(bin_levels, function(b) sum(sp_bins == b))

n_iter <- 10000
rand_rho <- numeric(n_iter)
for (i in seq_len(n_iter)) {
  used <- logical(length(genes))
  gs <- integer(0)
  for (k in seq_along(bin_levels)) {
    cand <- bin_pool[[k]]
    avail <- cand[!used[cand]]
    need <- n_need[k]
    if (length(avail) >= need) {
      pick <- sample(avail, need)
    } else {
      pick <- avail
      cand2 <- dec_pool[[bin_dec[k]]]
      avail2 <- cand2[!used[cand2] & !(cand2 %in% pick)]
      pick <- c(pick, sample(avail2, need - length(avail)))
    }
    used[pick] <- TRUE
    gs <- c(gs, pick)
  }
  rs <- colMeans(resid_m[gs, , drop = FALSE], na.rm = TRUE)
  rand_rho[i] <- cor(w, rs, method = "spearman")
}
rho_obs <- cor(w, splice_resid_score, method = "spearman")
p_emp <- (sum(abs(rand_rho) >= abs(rho_obs)) + 1) / (n_iter + 1)
p_emp_txt <- if (p_emp < 1e-4) "p < 1e-04" else paste("p =", signif(p_emp, 2))
saveRDS(rand_rho, "results/rna_protein_residual_null.rds")

# fraction of spliceosome proteins with positive median residual in LHS vs LLS
lhs_mean <- rowMeans(resid_m[sp, grp2 == "LHS", drop = FALSE], na.rm = TRUE)
lls_mean <- rowMeans(resid_m[sp, grp2 == "LLS", drop = FALSE], na.rm = TRUE)
delta <- lhs_mean - lls_mean

sink("results/rna_protein_residuals.txt", split = TRUE)
cat("=== Direct RNA-protein residual test (CPTAC, n =", length(samples), "patients) ===\n\n")
cat("Per-protein model (protein ~ mRNA + WEBS):\n")
cat("  spliceosome proteins with WEBS beta > 0, FDR<0.05:",
    sum(beta_df$is_spliceosome & beta_df$beta_webs > 0 & beta_df$q < 0.05, na.rm = TRUE),
    "/", sum(beta_df$is_spliceosome), "\n")
cat("  genome-wide:", sum(beta_df$beta_webs > 0 & beta_df$q < 0.05, na.rm = TRUE), "/", nrow(beta_df), "\n")
ft <- fisher.test(matrix(c(sum(beta_df$is_spliceosome & beta_df$beta_webs > 0 & beta_df$q < 0.05, na.rm = TRUE),
                           sum(beta_df$is_spliceosome & !(beta_df$beta_webs > 0 & beta_df$q < 0.05), na.rm = TRUE),
                           sum(!beta_df$is_spliceosome & beta_df$beta_webs > 0 & beta_df$q < 0.05, na.rm = TRUE),
                           sum(!beta_df$is_spliceosome & !(beta_df$beta_webs > 0 & beta_df$q < 0.05), na.rm = TRUE)), 2, byrow = TRUE))
cat("  Fisher enrichment of positive WEBS-beta in spliceosome: OR =", round(ft$estimate, 2),
    " p =", signif(ft$p.value, 3), "\n\n")
cat("Model B (complex residual score):\n")
cat("  Spliceosome residual score vs WEBS: rho =", round(ct$estimate, 3),
    " p =", signif(ct$p.value, 3), "\n")
cat("  ANOVA by group p =", signif(anova_p, 3), " KW p =", signif(kw$p.value, 3), "\n")
cat("  LHS-LLS residual delta: positive in", sum(delta > 0), "/", length(delta), "proteins\n")
cat("  Abundance+variance-matched random sets (10,000, non-spliceosome pool,\n")
cat("    without replacement): empirical", p_emp_txt, "\n")
sink()
cat(readLines("results/rna_protein_residuals.txt"), sep = "\n")

# figures
p1 <- ggplot(df, aes(group, resid_score, fill = group)) +
  geom_boxplot(outlier.size = 0.8) +
  scale_fill_manual(values = c("LLS" = "#1f77b4", "LDTS" = "#ff7f0e", "LHS" = "#d62728")) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(x = NULL, y = "Spliceosome residual score (protein − mRNA prediction)",
       title = "Spliceosome residual score by WEBS group",
       subtitle = paste0("KW p = ", signif(kw$p.value, 2), "; vs WEBS rho = ",
                         round(ct$estimate, 2), "; matched permutation ", p_emp_txt))
ggsave("figures_final/fig_rna_protein_residuals.pdf", p1, width = 6.5, height = 5)

beta_df$set <- ifelse(beta_df$is_spliceosome, "Spliceosome", "Non-spliceosome")
p2 <- ggplot(beta_df, aes(set, beta_webs, fill = set)) +
  geom_violin(scale = "width") + geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = c("Non-spliceosome" = "grey75", "Spliceosome" = "#7B4EA3")) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "WEBS beta (mRNA-adjusted)",
       title = "Per-protein WEBS effect after mRNA adjustment")
ggsave("figures_final/fig_residual_beta_dist.pdf", p2, width = 6, height = 5)

# combined two-panel supplementary figure (S36)
suppressMessages(library(patchwork))
fig_s36 <- p1 / p2 + plot_annotation(tag_levels = "a")
ggsave("figures_final/fig_S36_residuals.pdf", fig_s36, width = 6.5, height = 10)
ggsave("figures_final/fig_S36_residuals.png", fig_s36, width = 6.5, height = 10, dpi = 300)
cat("Done.\n")
