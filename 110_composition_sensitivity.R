# =============================================================================
# 110. Task 6 (V4.0): cell-composition sensitivity for KNCC module changes
#      Composition proxies (Neftel states + cell-type markers) added to the
#      competing-hypothesis model; overlap audit; classification of effects
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

b <- readRDS("results/kncc/dual_omics_bundle.rds")
delta <- b$prot_delta

# frozen modules
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f <- strsplit(l, "\t")[[1]]; f[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
hm <- lapply(hm, function(g) intersect(g, rownames(delta)))
mods <- list(
  PROLIF  = intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                        "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(delta)),
  NEURO   = intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                        "GAP43","TUBB3","MAP2"), rownames(delta)),
  ROS     = hm$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,
  IL2STAT5 = hm$HALLMARK_IL2_STAT5_SIGNALING,
  HYPOXIA = hm$HALLMARK_HYPOXIA,
  EMT     = hm$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  OXPHOS  = hm$HALLMARK_OXIDATIVE_PHOSPHORYLATION
)

# composition / state proxies (literature + Neftel 2019)
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
proxies <- lapply(proxies, function(g) intersect(g, rownames(delta)))
proxies <- proxies[sapply(proxies, length) >= 4]

sink("results/kncc/composition_sensitivity.txt", split = TRUE)
cat("=== Task 6: composition sensitivity ===\n\n")
cat("--- overlap audit (proxy genes vs modules) ---\n")
for (pn in names(proxies)) {
  ov <- sapply(mods, function(m) length(intersect(proxies[[pn]], m)))
  if (any(ov > 0)) cat(pn, "overlaps:", paste(names(ov)[ov > 0], ov[ov > 0], collapse = "; "), "\n")
}

# per-patient proxy deltas
px_delta <- t(sapply(proxies, function(g) colMeans(delta[g, , drop = FALSE], na.rm = TRUE)))
colnames(px_delta) <- colnames(delta)

# adjudication with composition predictors
d_pro <- colMeans(delta[setdiff(mods$PROLIF, unlist(proxies)), , drop = FALSE], na.rm = TRUE)
d_neu <- colMeans(delta[setdiff(mods$NEURO, unlist(proxies)), , drop = FALSE], na.rm = TRUE)
d_imm <- colMeans(delta[intersect(c("CD14","AIF1","FCGR3A","CD68","CD163","ITGAM","PTPRC",
                                    "CSF1R","TMEM119","CX3CR1"), rownames(delta)), , drop = FALSE], na.rm = TRUE)

cat("\n--- module adjudication: without vs with composition proxies ---\n")
res <- lapply(names(mods), function(mn) {
  g <- mods[[mn]]
  d_m <- colMeans(delta[g, , drop = FALSE], na.rm = TRUE)
  ok <- !is.na(d_m)
  # base model (as frozen)
  f1 <- lm(d_m[ok] ~ d_pro[ok] + d_neu[ok] + d_imm[ok])
  b1 <- summary(f1)$coefficients[1, ]
  # extended model with composition proxies (excluding genes of target module)
  px_use <- lapply(proxies, function(pg) setdiff(pg, g))
  px_use <- px_use[sapply(px_use, length) >= 4]
  X <- t(sapply(px_use, function(pg) colMeans(delta[pg, , drop = FALSE], na.rm = TRUE)))
  df2 <- data.frame(y = d_m[ok], p = d_pro[ok], n = d_neu[ok], i = d_imm[ok], t(X[, ok, drop = FALSE]))
  f2 <- lm(y ~ ., data = df2)
  b2 <- summary(f2)$coefficients[1, ]
  data.frame(module = mn,
             base_intercept = b1[1], base_p = b1[4],
             comp_intercept = b2[1], comp_p = b2[4],
             base_r2 = summary(f1)$r.squared, comp_r2 = summary(f2)$r.squared)
})
res <- do.call(rbind, res)
res$class <- ifelse(res$comp_p < 0.05, "stable residual (composition-adjusted)",
             ifelse(res$base_p < 0.05, "composition-associated (lost after comp. adjustment)",
                    "state-coupled / n.s."))
print(res, row.names = FALSE)
sink()
write.csv(res, "results/kncc/composition_sensitivity.csv", row.names = FALSE)

# figure
res$module <- factor(res$module, levels = res$module[order(res$comp_intercept)])
p <- ggplot(res, aes(module, comp_intercept)) +
  geom_col(aes(fill = comp_p < 0.05)) +
  geom_point(aes(y = base_intercept), color = "grey30", shape = 4, size = 3) +
  coord_flip() + geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("TRUE" = "#7B4EA3", "FALSE" = "grey75"),
                    labels = c("n.s.", "p<0.05")) +
  theme_bw(base_size = 11) +
  labs(y = "Intercept (state-independent module change)", x = NULL, fill = "composition-adj.",
       title = "Module effects before (x) and after (bar) composition adjustment",
       caption = "bar = adjusted for state + cell-composition proxies; x = state-only model")
ggsave("figures_kncc/composition_sensitivity.pdf", p, width = 7.5, height = 5)
cat("Done.\n")
