# =============================================================================
# 42. WGCNA on CPTAC proteome + acetylome: module-WEBS association
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(WGCNA); library(dplyr); library(ggplot2) })
options(stringsAsFactors = FALSE)
set.seed(42)

webs <- read.csv("results/cptac_webs.csv")
prot <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct",
                             row.names = 1, check.names = FALSE))
prot <- prot[rowSums(!is.na(prot)) > ncol(prot) * 0.6, ]
prot <- prot[order(apply(prot, 1, var, na.rm = TRUE), decreasing = TRUE)[1:4000], ]

common <- intersect(colnames(prot), webs$sample)
prot <- prot[, common]
w <- webs$WEBS[match(common, webs$sample)]

cat("WGCNA on", nrow(prot), "proteins x", ncol(prot), "samples\n")

datExpr <- t(prot)
sft <- pickSoftThreshold(datExpr, verbose = 0)
power <- sft$powerEstimate
if (is.na(power)) power <- 6
cat("Soft power:", power, "\n")

net <- blockwiseModules(datExpr, power = power, TOMType = "unsigned",
                        minModuleSize = 30, reassignThreshold = 0,
                        mergeCutHeight = 0.25, numericLabels = TRUE,
                        pamRespectsDendro = FALSE, verbose = 0)

MEs <- moduleEigengenes(datExpr, colors = net$colors)$eigengenes
mod_trait <- cor(MEs, w, use = "p")
mod_p <- corPvalueStudent(mod_trait, nrow(datExpr))

mt_df <- data.frame(module = colnames(MEs), correlation = as.vector(mod_trait),
                    p = as.vector(mod_p), q = p.adjust(as.vector(mod_p), method = "BH"))
mt_df <- mt_df[order(-abs(mt_df$correlation)), ]
write.csv(mt_df, "results/wgcna_module_webs.csv", row.names = FALSE)

sink("results/wgcna_summary.txt", split = TRUE)
cat("=== WGCNA module-WEBS association ===\n")
cat("Modules:", length(unique(net$colors)), "\n\n")
print(mt_df)
sink()

# top module: gene list + enrichment marker
top_mod <- mt_df$module[1]
top_genes <- colnames(datExpr)[net$colors == as.integer(sub("ME", "", top_mod))]
writeLines(top_genes, "results/wgcna_top_module_genes.txt")

# module-trait heatmap
library(pheatmap)
pdf("figures_final/fig_wgcna_moduletrait.pdf", width = 7, height = 5)
mat <- t(mod_trait); rownames(mat) <- "WEBS"
colnames(mat) <- colnames(MEs)
pheatmap(mat, display_numbers = matrix(sprintf("%.2f", mod_trait), nrow = 1),
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         main = "Module eigengene - WEBS correlation")
dev.off()
cat("Done. Top module:", top_mod, "with", length(top_genes), "genes, r =",
    round(mt_df$correlation[1], 3), "q =", signif(mt_df$q[1], 3), "\n")
