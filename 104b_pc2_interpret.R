setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(dplyr))
b <- readRDS("results/kncc/dual_omics_bundle.rds")
delta <- b$prot_delta
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f <- strsplit(l, "\t")[[1]]; f[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
keep_mods <- c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_DNA_REPAIR",
               "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS",
               "HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
               "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_INTERFERON_GAMMA_RESPONSE",
               "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
               "HALLMARK_ANGIOGENESIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_PI3K_AKT_MTOR_SIGNALING",
               "HALLMARK_P53_PATHWAY","HALLMARK_APOPTOSIS","HALLMARK_HYPOXIA",
               "HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_PROTEIN_SECRETION",
               "HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_MYC_TARGETS_V1")
keep_mods <- keep_mods[keep_mods %in% names(hm)]
mods <- lapply(hm[keep_mods], function(g) intersect(g, rownames(delta)))
mods$NEURONAL_CURATED <- intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2","GAP43","TUBB3","MAP2"), rownames(delta))
mods$PROLIF_CURATED <- intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(delta))
mods <- mods[sapply(mods, length) >= 8]
Md <- t(sapply(mods, function(g) colMeans(delta[g, , drop = FALSE], na.rm = TRUE)))
Md <- Md[, colSums(is.na(Md)) == 0]
pc <- prcomp(t(Md), scale. = TRUE)
ld <- round(pc$rotation[, 1:3], 3)
ld <- ld[order(-abs(ld[, 2])), ]
cat("PC2 loadings (sorted by |loading|):\n")
print(ld[, 1:3])
# correlations of PC2 with each module delta
cat("\nPC2 vs module deltas (Spearman):\n")
cr <- apply(Md, 1, function(v) cor(pc$x[, 2], v, method = "spearman"))
print(round(sort(cr), 3))
