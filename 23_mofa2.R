# =============================================================================
# 23. MOFA2 multi-omics integration on CPTAC-GBM (4 views, 99 samples)
# =============================================================================
# Views: proteome, acetylome (site->gene mean), phosphoproteome (site->gene
#        mean), RNAseq — same 99 samples
# Goal: find factor(s) tracking WEBS; check spliceosome loadings in protein view
# Output: results/mofa2_summary.txt, results/mofa2_factors.csv
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(MOFA2); library(dplyr); library(ggplot2) })
set.seed(42)

# point reticulate to the mofa conda env (mofapy2 backend)
reticulate::use_python("E:/Miniconda3/envs/mofa/python.exe", required = TRUE)

webs <- read.csv("results/cptac_webs.csv", stringsAsFactors = FALSE)

load_view <- function(f, collapse = FALSE) {
  m <- as.matrix(read.delim(f, row.names = 1, check.names = FALSE))
  m <- m[rowSums(!is.na(m)) > ncol(m) * 0.5, ]  # keep features quantified in >50%
  if (collapse) {
    g <- sub("[-_\\.][A-Za-z]*[0-9]+.*$", "", rownames(m))  # strip site suffix
    m <- rowsum(m, g, reorder = FALSE, na.rm = TRUE) /
         pmax(1, as.vector(table(g)[rownames(rowsum(m, g, reorder = FALSE))]))
  }
  m
}

prot <- load_view("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct")
acet <- load_view("data/raw/CPTAC/HS_CPTAC_GBM_acetylome_mssm_site.cct", collapse = TRUE)
phos <- load_view("data/raw/CPTAC/HS_CPTAC_GBM_phosphoproteome_mssm_site.cct", collapse = TRUE)
rna  <- load_view("data/raw/CPTAC/HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct")

# top variable features per view (MOFA speed)
topv <- function(m, n = 3000) m[order(apply(m, 1, var, na.rm = TRUE), decreasing = TRUE)[1:min(n, nrow(m))], ]

views <- list(proteome = topv(prot), acetylome = topv(acet),
              phosphoproteome = topv(phos), rnaseq = topv(rna))
samples <- Reduce(intersect, lapply(views, colnames))
views <- lapply(views, function(m) m[, samples])
cat("Views x samples:", paste(sapply(views, nrow), collapse = "/"), "x", length(samples), "\n")

mofa <- create_mofa(views)
opts <- get_default_model_options(mofa); opts$num_factors <- 10
train <- get_default_training_options(mofa); train$seed <- 42
mofa <- prepare_mofa(mofa, model_options = opts, training_options = train)
mofa <- run_mofa(mofa, use_basilisk = FALSE)

factors <- get_factors(mofa)$group1
colnames(factors) <- paste0("Factor", seq_len(ncol(factors)))
w <- webs$WEBS[match(samples, webs$sample)]

fres <- t(sapply(seq_len(ncol(factors)), function(i) {
  ct <- cor.test(factors[, i], w, method = "spearman")
  c(factor = i, rho = unname(ct$estimate), p = ct$p.value)
}))
fres <- as.data.frame(fres); fres$q <- p.adjust(fres$p, method = "BH")
write.csv(fres, "results/mofa2_factor_webs_corr.csv", row.names = FALSE)

# top factor: spliceosome loadings in protein view
top_f <- fres$factor[which.max(abs(fres$rho))]
lw <- get_weights(mofa, views = "proteome", factors = top_f)[[1]][, 1]
source("08_cptac_proteomics.R")
sp <- get_spliceosome_genes()
in_sp <- names(lw) %in% sp
wt <- wilcox.test(lw[in_sp], lw[!in_sp])

var_exp <- get_variance_explained(mofa)$r2_per_factor$group1

sink("results/mofa2_summary.txt")
cat("=== MOFA2 CPTAC 4-view integration ===\n\n")
cat("Samples:", length(samples), "| Factors: 10\n\n")
cat("Variance explained per factor (%):\n")
print(round(var_exp, 2))
cat("\nFactor-WEBS correlations (Spearman):\n")
print(fres)
cat("\nTop WEBS-associated factor: Factor", top_f, " (rho =",
    round(fres$rho[fres$factor == top_f], 3), ")\n")
cat("Spliceosome protein loadings vs others (Wilcoxon): p =", signif(wt$p.value, 3), "\n")
cat("Median loading spliceosome:", round(median(lw[in_sp]), 4),
    " vs others:", round(median(lw[!in_sp]), 4), "\n")
sink()
cat(readLines("results/mofa2_summary.txt"), sep = "\n")

factors_df <- data.frame(sample = samples, factors, WEBS = w)
write.csv(factors_df, "results/mofa2_factors.csv", row.names = FALSE)
saveRDS(mofa, "results/mofa2_model.rds")
message("Done.")
