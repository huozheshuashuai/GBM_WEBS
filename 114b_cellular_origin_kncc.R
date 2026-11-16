# =============================================================================
# 114b. Cellular origin of the second axis, KNCC-proteome-based
#       Correlate T1 module scores (ROS/IL2STAT5/EMT + refs) with cell-type
#       proxy protein scores across 109 patients (cross-sectional)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

b <- readRDS("results/kncc/dual_omics_bundle.rds")
pm <- b$prot_samples; cli <- b$prot_annot
t1 <- cli$sample[cli$timepoint == "T1"]
pm1 <- pm[, t1]

gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f2 <- strsplit(l, "\t")[[1]]; f2[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
hm <- lapply(hm, function(g) intersect(g, rownames(pm1)))
mods <- list(
  PROLIF  = intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                        "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(pm1)),
  NEURO   = intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                        "GAP43","TUBB3","MAP2"), rownames(pm1)),
  ROS     = hm$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,
  IL2STAT5 = hm$HALLMARK_IL2_STAT5_SIGNALING,
  HYPOXIA = hm$HALLMARK_HYPOXIA,
  EMT     = hm$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  OXPHOS  = hm$HALLMARK_OXIDATIVE_PHOSPHORYLATION
)
proxies <- list(
  MES_state   = c("CD44","VIM","S100A11","S100A10","ANXA1","LGALS1","MT2A","C1S","NAMPT","EMP3"),
  AC_state    = c("GFAP","S100B","ALDOC","SLC1A3","AQP4","CLU","MLC1","ATP1A2","FABP7","AGT"),
  OPC_state   = c("SOX10","PLP1","BCAS1","GPR17","MYRF","MAG","MOG","CNP","MOBP","OPALIN"),
  NPC_state   = c("SOX4","SOX11","STMN1","STMN2","DCX","DLL3","BEX1","GAP43","SNAP25","MARCKSL1"),
  Microglia   = c("PTPRC","CD14","CSF1R","C1QA","AIF1","TMEM119","P2RY12","CX3CR1","CD68","CD163"),
  Endothelial = c("PECAM1","VWF","CLDN5","CDH5","KDR","ESAM","FLT1","ENG"),
  Pericyte    = c("PDGFRB","RGS5","ACTA2","CSPG4","MCAM","DES"),
  Oligodendrocyte = c("MBP","PLP1","MOG","MAG","OPALIN","MOBP","CNP","TF"),
  Neuron      = c("RBFOX3","SNAP25","SYP","SYT1","STMN2","TUBB3","MAP2","NEFL")
)
proxies <- lapply(proxies, function(g) intersect(g, rownames(pm1)))
proxies <- proxies[sapply(proxies, length) >= 4]

S <- t(sapply(mods, function(g) colMeans(pm1[g, , drop = FALSE], na.rm = TRUE)))
P <- t(sapply(proxies, function(g) colMeans(pm1[g, , drop = FALSE], na.rm = TRUE)))

res <- expand.grid(module = rownames(S), proxy = rownames(P))
res$rho <- mapply(function(a, b) cor(S[a, ], P[b, ], method = "spearman", use = "complete.obs"),
                  res$module, res$proxy)
res$p <- mapply(function(a, b) suppressWarnings(cor.test(S[a, ], P[b, ], method = "spearman")$p.value),
                res$module, res$proxy)
res$q <- p.adjust(res$p, "BH")

sink("results/kncc/cellular_origin_kncc.txt", split = TRUE)
cat("=== cellular origin: second-axis modules vs cell-type proxies (KNCC T1, n =", ncol(pm1), ") ===\n\n")
for (mn in rownames(S)) {
  d <- res[res$module == mn, ]
  d <- d[order(-d$rho), ]
  cat(mn, ": top proxies:", paste0(d$proxy[1:3], " (", round(d$rho[1:3], 2), ")", collapse = ", "), "\n")
}
cat("\n")
sig <- res[res$q < 0.05, ]
print(sig[order(sig$module, -sig$rho), ], row.names = FALSE)
sink()

res$module <- factor(res$module, levels = rownames(S))
res$proxy <- factor(res$proxy, levels = rownames(P))
p <- ggplot(res, aes(proxy, module, fill = rho)) +
  geom_tile() +
  geom_text(aes(label = ifelse(q < 0.05, "*", "")), size = 5) +
  scale_fill_gradient2(low = "#3D9B8C", mid = "white", high = "#7B4EA3", midpoint = 0) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "cell-type proxy (protein markers)", y = NULL, fill = "Spearman rho",
       title = "Second-axis modules vs cell-type proxies (KNCC primary tumors)",
       caption = "* FDR < 0.05")
ggsave("figures_kncc/cellular_origin_kncc.pdf", p, width = 8, height = 5)
write.csv(res, "results/kncc/cellular_origin_kncc.csv", row.names = FALSE)
cat("Done.\n")
