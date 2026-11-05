# =============================================================================
# 62. Threshold sensitivity of the 107/111 spliceosome result
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(limma); library(dplyr) })
source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()

prot <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct", check.names = FALSE)
pmat <- as.matrix(sapply(prot[, -1], as.numeric)); rownames(pmat) <- prot[[1]]
webs <- read.csv("results/cptac_webs.csv")
common <- intersect(colnames(pmat), webs$sample)
pmat <- pmat[, common]
w <- webs$WEBS[match(common, webs$sample)]

sink("results/threshold_sensitivity.txt", split = TRUE)
cat("=== Threshold sensitivity: spliceosome up/down in high vs low WEBS ===\n\n")
settings <- list(c(0.20, 0.60), c(0.25, 0.50), c(0.30, 0.40), c(0.50, 0.00))
for (s in settings) {
  lo_q <- s[1]; hi_q <- 1 - s[1]
  if (s[2] == 0) {  # median split
    grp <- ifelse(w > median(w), "high", "low")
  } else {
    grp <- ifelse(w <= quantile(w, lo_q), "low",
                  ifelse(w >= quantile(w, hi_q), "high", NA))
  }
  sel <- !is.na(grp)
  g <- factor(grp[sel], levels = c("low", "high"))
  des <- model.matrix(~ g)
  sp <- intersect(splice, rownames(pmat))
  fit <- eBayes(lmFit(pmat[sp, sel], des))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  n_up <- sum(tt$logFC > 0); n_up_sig <- sum(tt$logFC > 0 & tt$adj.P.Val < 0.05)
  label <- ifelse(s[2] == 0, "median split",
                  sprintf("bottom/top %d%%", round(lo_q * 100)))
  cat(sprintf("%-18s n_low=%2d n_high=%2d  spliceosome up: %d/%d (FDR<0.05: %d)  down: %d\n",
              label, sum(g == "low"), sum(g == "high"), n_up, nrow(tt), n_up_sig, sum(tt$logFC < 0)))
}
sink()
cat(readLines("results/threshold_sensitivity.txt"), sep = "\n")
