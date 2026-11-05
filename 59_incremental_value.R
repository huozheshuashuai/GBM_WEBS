# =============================================================================
# 59. Incremental value of WEBS beyond established GBM cell-state programs
#     (i) CPTAC protein level: spliceosome ~ WEBS | proliferation, MES, hypoxia
#     (ii) scRNA: WEBS vs spliceosome expression within each Neftel state
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(Seurat); library(ppcor); library(ggplot2); library(dplyr) })
source("08_cptac_proteomics.R")
source("tools/WEBS/webs.R")

proteome <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct",
                                 row.names = 1, check.names = FALSE))
webs <- read.csv("results/cptac_webs.csv")
cs <- intersect(colnames(proteome), webs$sample)
pw <- webs$WEBS[match(cs, webs$sample)]

splice_in <- intersect(get_spliceosome_genes(), rownames(proteome))
score_p <- colMeans(proteome[splice_in, cs, drop = FALSE], na.rm = TRUE)

mk <- function(genes) colMeans(proteome[intersect(genes, rownames(proteome)), cs, drop = FALSE], na.rm = TRUE)
prolif <- mk(c("MKI67","TOP2A","PCNA","MCM2","MCM6"))
mes    <- mk(c("VIM","CD44","FN1","SERPINE1","TGFBI"))
hypox  <- mk(c("CA9","VEGFA","LDHA","SLC2A1"))
opc    <- mk(c("OLIG2","SOX10","PDGFRA"))

sink("results/incremental_value.txt", split = TRUE)
cat("=== WEBS incremental value (CPTAC protein level) ===\n\n")
cat(sprintf("Spliceosome ~ WEBS:                      rho = %.3f\n",
            cor(pw, score_p, method = "spearman")))
for (nm in c("prolif","mes","hypox","opc")) {
  v <- get(nm)
  pc <- pcor.test(pw, score_p, v, method = "spearman")
  cat(sprintf("  controlling %-8s                 rho = %.3f, p = %.3g\n", nm, pc$estimate, pc$p.value))
}
pc_all <- pcor.test(pw, score_p, cbind(prolif, mes, hypox, opc), method = "spearman")
cat(sprintf("  controlling ALL four                   rho = %.3f, p = %.3g\n", pc_all$estimate, pc_all$p.value))

# WEBS vs each state score (is WEBS just a state proxy?)
cat("\nWEBS correlation with state scores (protein level):\n")
for (nm in c("prolif","mes","hypox","opc")) {
  ct <- cor.test(pw, get(nm), method = "spearman")
  cat(sprintf("  WEBS vs %-8s rho = %+.3f, p = %.3g\n", nm, ct$estimate, ct$p.value))
}

# --- scRNA: within-state association ---
mal <- readRDS("data/processed/GSE138794_malignant.rds")
md <- mal@meta.data
splice_sc <- intersect(get_spliceosome_genes(), rownames(mal))
splice_score <- colMeans(GetAssayData(mal, layer = "data")[splice_sc, ])
md$splice_score <- as.numeric(splice_score[rownames(md)])
md$Neftel_state <- apply(md[, c("Neftel_MES","Neftel_AC","Neftel_NPC","Neftel_OPC")], 1,
                         function(x) c("MES","AC","NPC","OPC")[which.max(x)])
cat("\nscRNA: WEBS vs spliceosome expression within each Neftel state (per-patient median rho):\n")
for (st in c("MES","AC","NPC","OPC")) {
  sub <- md[md$Neftel_state == st, ]
  per_pat <- sapply(split(sub, sub$sample), function(d) {
    if (nrow(d) < 100 || sd(d$WEBS) == 0) return(NA)
    cor(d$WEBS, d$splice_score, method = "spearman")
  })
  per_pat <- per_pat[!is.na(per_pat)]
  if (length(per_pat) >= 3) {
    wt <- wilcox.test(per_pat, mu = 0)
    cat(sprintf("  %s: median rho = %+.3f, n patients = %d, Wilcoxon p = %.3g\n",
                st, median(per_pat), length(per_pat), wt$p.value))
  } else {
    cat(sprintf("  %s: only %d patients with enough cells\n", st, length(per_pat)))
  }
}
sink()
cat(readLines("results/incremental_value.txt"), sep = "\n")
cat("Done.\n")
