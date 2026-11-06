# =============================================================================
# 64. CPTAC six-gene survival + GLASS first-pair sensitivity
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(survival); library(dplyr) })

PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)

sink("results/cptac_glass_sensitivity.txt", split = TRUE)

# --- CPTAC signature survival ---
cat("=== CPTAC-GBM: six-gene signature survival ===\n")
rna <- read.delim("E:/KIMI/data/HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct", check.names = FALSE)
m <- as.matrix(sapply(rna[, -1], as.numeric)); rownames(m) <- rna[[1]]
m <- m[!is.na(rownames(m)), ]
g <- intersect(names(PAPER_COEFS), rownames(m))
risk <- as.vector(t(m[g, ]) %*% PAPER_COEFS[g])
cli <- read.delim("E:/KIMI/data/HS_CPTAC_GBM_CLI.tsi", check.names = FALSE, stringsAsFactors = FALSE)
cat("CLI cols:", paste(head(colnames(cli), 15), collapse = ", "), "\n")
cat("CLI rows:", nrow(cli), "\n")
# find survival fields
print(head(cli, 3))
sink()
