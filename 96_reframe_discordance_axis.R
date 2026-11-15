# =============================================================================
# 96. WEBS-free rebuild: CPTAC spliceosome protein-mRNA discordance axis
#     Per-patient spliceosome protein score vs mRNA score;
#     per-protein residuals; discordance score = mean residual (protein~mRNA);
#     abundance+variance-matched permutation control (strict, from fixed 80)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2); library(patchwork) })
set.seed(42)
dir.create("results/reframe", showWarnings = FALSE)
dir.create("figures_reframe", showWarnings = FALSE)

prot <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct", check.names = FALSE)
rna  <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct", check.names = FALSE)
pm <- as.matrix(sapply(prot[, -1], as.numeric)); rownames(pm) <- prot[[1]]
rm_ <- as.matrix(sapply(rna[, -1], as.numeric)); rownames(rm_) <- rna[[1]]

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()

genes <- intersect(rownames(pm), rownames(rm_))
samples <- intersect(colnames(pm), colnames(rm_))
pm <- pm[genes, samples]; rm_ <- rm_[genes, samples]
sp <- intersect(splice, genes)
cat("CPTAC paired samples:", length(samples), "| spliceosome proteins:", length(sp), "\n")

# ---- per-patient complex-level discordance ----
prot_score <- colMeans(pm[sp, ], na.rm = TRUE)
mrna_score <- colMeans(rm_[sp, ], na.rm = TRUE)
ct <- cor.test(prot_score, mrna_score, method = "spearman")
cat("spliceosome protein score vs mRNA score: rho =", round(ct$estimate, 3),
    " p =", signif(ct$p.value, 3), "\n")

# ---- per-protein residuals (protein ~ mRNA across patients) ----
resid_m <- matrix(NA_real_, length(genes), length(samples),
                  dimnames = list(genes, samples))
for (g in genes) {
  fit <- tryCatch(lm(pm[g, ] ~ rm_[g, ], na.action = na.exclude), error = function(e) NULL)
  if (!is.null(fit)) resid_m[g, ] <- resid(fit)
}
spds <- colMeans(resid_m[sp, ], na.rm = TRUE)   # spliceosome protein discordance score per patient
cat("SPDS range:", round(range(spds), 2), "\n")

# per-protein: fraction with positive residual association to SPDS
get_beta <- function(g) {
  fit <- tryCatch(lm(pm[g, ] ~ rm_[g, ] + spds), error = function(e) NULL)
  if (is.null(fit)) return(c(NA, NA))
  co <- summary(fit)$coefficients
  if (nrow(co) < 3) return(c(NA, NA))
  c(co["spds", "Estimate"], co["spds", "Pr(>|t|)"])
}
ba <- t(sapply(genes, get_beta))
beta_df <- data.frame(gene = genes, beta = ba[, 1], p = ba[, 2])
beta_df$q <- p.adjust(beta_df$p, "BH")
beta_df$is_spliceosome <- beta_df$gene %in% sp
write.csv(beta_df, "results/reframe/cptac_spds_perprotein.csv", row.names = FALSE)

n_pos <- sum(beta_df$is_spliceosome & beta_df$beta > 0 & beta_df$q < 0.05, na.rm = TRUE)
n_sp_tested <- sum(beta_df$is_spliceosome)
cat("spliceosome proteins positively tracking SPDS (FDR<0.05):", n_pos, "/", n_sp_tested, "\n")
cat("(expected tautology check: SPDS is built from these proteins' residuals;\n")
cat(" the NON-circular result is what tracks SPDS genome-wide)\n")
cat("genome-wide positive:", sum(beta_df$beta > 0 & beta_df$q < 0.05, na.rm = TRUE), "\n")

# genome-wide: which non-spliceosome proteins track the discordance axis?
top_nonsp <- beta_df %>% filter(!is_spliceosome, q < 0.05) %>%
  arrange(desc(abs(beta))) %>% head(20)
print(top_nonsp[, c("gene", "beta", "q")], row.names = FALSE)

# ---- strict permutation: is the SPDS spread (inter-patient variance of residual score)
# larger than matched random complexes? ----
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
pool <- setdiff(genes, sp)
sp_bins <- bin[sp]; bin_levels <- sort(unique(sp_bins))
bin_pool <- lapply(bin_levels, function(b) which(bin == b & names(bin) %in% pool))
names(bin_pool) <- bin_levels
dec_pool <- lapply(1:10, function(d) which(dec == d & names(dec) %in% pool))
bin_dec <- ceiling(bin_levels / 3)
n_need <- sapply(bin_levels, function(b) sum(sp_bins == b))

obs_var <- var(spds)
n_iter <- 10000
rand_var <- numeric(n_iter)
for (i in seq_len(n_iter)) {
  used <- logical(length(genes)); gs <- integer(0)
  for (k in seq_along(bin_levels)) {
    cand <- bin_pool[[k]]; avail <- cand[!used[cand]]; need <- n_need[k]
    if (length(avail) >= need) {
      pick <- sample(avail, need)
    } else {
      pick <- avail
      cand2 <- dec_pool[[bin_dec[k]]]
      avail2 <- cand2[!used[cand2] & !(cand2 %in% pick)]
      pick <- c(pick, sample(avail2, need - length(avail)))
    }
    used[pick] <- TRUE; gs <- c(gs, pick)
  }
  rand_var[i] <- var(colMeans(resid_m[gs, , drop = FALSE], na.rm = TRUE))
}
p_emp <- (sum(rand_var >= obs_var) + 1) / (n_iter + 1)
cat("\nSPDS inter-patient variance:", round(obs_var, 4),
    "| matched random sets (10,000): empirical p =", signif(p_emp, 3), "\n")

saveRDS(list(spds = spds, prot_score = prot_score, mrna_score = mrna_score,
             resid_m = resid_m[sp, ], samples = samples, beta_df = beta_df),
        "results/reframe/cptac_discordance_axis.rds")

# ---- figures ----
df <- data.frame(prot = prot_score, mrna = mrna_score, spds = spds)
p1 <- ggplot(df, aes(mrna, prot)) +
  geom_point(color = "#7B4EA3", alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  theme_bw(base_size = 12) +
  labs(x = "Spliceosome mRNA score", y = "Spliceosome protein score",
       title = sprintf("CPTAC-GBM (n = %d): rho = %.2f, p = %s",
                       length(samples), ct$estimate, signif(ct$p.value, 2)),
       subtitle = "protein and transcript layers are decoupled at complex level")
p2 <- ggplot(data.frame(x = rand_var), aes(x)) +
  geom_histogram(bins = 50, fill = "grey75") +
  geom_vline(xintercept = obs_var, color = "#7B4EA3", linewidth = 1) +
  theme_bw(base_size = 12) +
  labs(x = "Inter-patient variance of residual score (matched random sets)",
       y = "Count",
       title = sprintf("Spliceosome discordance spread: var = %.3f, empirical p = %s",
                       obs_var, signif(p_emp, 2)))
ggsave("figures_reframe/cptac_discordance.pdf", p1 | p2, width = 11, height = 5)
cat("Done.\n")
