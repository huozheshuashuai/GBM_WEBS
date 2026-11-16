# =============================================================================
# 118. GeoMx spatial localization of the second axis (GSE254875)
#      79 DCC files (WTA) -> gene x ROI counts; frozen KNCC signatures;
#      CD64 (myeloid) vs GFAP (neural) ROIs, primary vs recurrent;
#      patient as statistical unit
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2); library(tidyr) })
set.seed(42)
dir.create("results/geomx", showWarnings = FALSE)

# ---------- RTS -> gene map ----------
rts_map <- read.delim("data/raw/wta_rts_gene_map.tsv")
rts2gene <- setNames(rts_map$Gene, rts_map$RTS_ID)

# ---------- parse all DCC files ----------
dcc_files <- list.files("E:/KIMI/data/GSE254875_dcc", pattern = "\\.dcc\\.gz$", full.names = TRUE)
cat("DCC files:", length(dcc_files), "\n")

parse_dcc <- function(f) {
  lines <- readLines(gzfile(f))
  cs <- grep("^<Code_Summary>", lines) + 1
  ce <- grep("^</Code_Summary>", lines) - 1
  if (!length(cs) || !length(ce) || ce < cs) return(NULL)
  tab <- lines[cs:ce]
  tab <- tab[grepl("^RTS", tab)]
  cnt <- as.numeric(sub(".*,", "", tab))
  names(cnt) <- sub(",.*", "", tab)
  cnt
}
counts_list <- lapply(dcc_files, parse_dcc)
ok <- !sapply(counts_list, is.null)
dcc_files <- dcc_files[ok]; counts_list <- counts_list[ok]
cat("parsed:", length(counts_list), "\n")

# align on the RTS union
all_rts <- unique(unlist(lapply(counts_list, names)))
cmat <- matrix(0, length(all_rts), length(counts_list),
               dimnames = list(all_rts, sub("\\.dcc\\.gz$", "", basename(dcc_files))))
for (i in seq_along(counts_list)) cmat[names(counts_list[[i]]), i] <- counts_list[[i]]
cat("raw RTS x ROI matrix:", nrow(cmat), "x", ncol(cmat), "\n")

# collapse to gene level
gene <- rts2gene[rownames(cmat)]
keep <- !is.na(gene)
cmat <- cmat[keep, ]; gene <- gene[keep]
gmat <- rowsum(cmat, gene)
cat("gene-level matrix:", nrow(gmat), "genes x", ncol(gmat), "ROIs\n")

# ---------- ROI metadata from GPL24676 series matrix ----------
sm <- readLines(gzfile("E:/KIMI/data/GSE254875-GPL24676_series_matrix.txt.gz"))
get_field <- function(tag) {
  l <- sm[grepl(paste0("^!", tag), sm)][1]
  strsplit(gsub('"', "", l), "\t")[[1]][-1]
}
gsm <- get_field("Sample_geo_accession")
title <- get_field("Sample_title")
roi_meta <- data.frame(gsm = gsm, title = title)
roi_meta$dcc_stub <- sapply(roi_meta$gsm, function(g) {
  hit <- colnames(cmat)[grepl(paste0("^", g, "_"), colnames(cmat))]
  if (length(hit)) hit[1] else NA_character_
})
roi_meta <- roi_meta[!is.na(roi_meta$dcc_stub), ]
cat("ROIs matched to DCC:", nrow(roi_meta), "\n")

# parse title: patient, timepoint, stain marker
roi_meta$patient <- sub("^([0-9]+)_.*", "\\1", roi_meta$title)
roi_meta$timepoint <- ifelse(grepl("primary", roi_meta$title), "Primary",
                      ifelse(grepl("recurrent", roi_meta$title), "Recurrent", NA))
roi_meta$stain <- ifelse(grepl("CD64", roi_meta$title), "CD64",
                  ifelse(grepl("GFAP", roi_meta$title), "GFAP", NA))
cat("\nROI x timepoint x stain:\n")
print(table(roi_meta$timepoint, roi_meta$stain))
cat("patients:", length(unique(roi_meta$patient)), "\n")
roi_meta <- roi_meta[!is.na(roi_meta$timepoint) & !is.na(roi_meta$stain), ]

# ---------- normalization: Q3 + log ----------
gmat <- gmat[, roi_meta$dcc_stub]
lib <- colSums(gmat)
q3 <- apply(gmat, 2, function(x) quantile(x[x > 0], 0.75, na.rm = TRUE))
norm <- t(t(gmat) / (q3 + 1) * 1000)
lnorm <- log2(norm + 0.5)

# ---------- frozen module signatures ----------
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f2 <- strsplit(l, "\t")[[1]]; f2[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
hm <- lapply(hm, function(g) intersect(g, rownames(lnorm)))
mods <- list(
  PROLIF  = intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                        "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(lnorm)),
  NEURO   = intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                        "GAP43","TUBB3","MAP2"), rownames(lnorm)),
  ROS     = hm$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,
  IL2STAT5 = hm$HALLMARK_IL2_STAT5_SIGNALING,
  HYPOXIA = hm$HALLMARK_HYPOXIA,
  EMT     = hm$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  OXPHOS  = hm$HALLMARK_OXIDATIVE_PHOSPHORYLATION
)
cat("\nmodule coverage in GeoMx:\n")
for (mn in names(mods)) cat(" ", mn, ":", length(mods[[mn]]), "\n")
proxies <- list(
  Microglia = intersect(c("PTPRC","CD14","CSF1R","C1QA","AIF1","TMEM119","P2RY12",
                          "CX3CR1","CD68","CD163"), rownames(lnorm)),
  Endothelial = intersect(c("PECAM1","VWF","CLDN5","CDH5","KDR","ESAM","FLT1","ENG"), rownames(lnorm)),
  MES_state = intersect(c("CD44","VIM","S100A11","S100A10","ANXA1","LGALS1","MT2A","C1S","NAMPT","EMP3"), rownames(lnorm))
)
cat("proxies:", paste(names(proxies), sapply(proxies, length), collapse = ", "), "\n")

# frozen PC2 axis score: use KNCC frozen loadings (from script 116 recompute)
b <- readRDS("results/kncc/dual_omics_bundle.rds")
kdelta <- b$prot_delta
keep_mods_pc <- c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_DNA_REPAIR",
                  "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS",
                  "HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                  "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_INTERFERON_GAMMA_RESPONSE",
                  "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                  "HALLMARK_ANGIOGENESIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_PI3K_AKT_MTOR_SIGNALING",
                  "HALLMARK_P53_PATHWAY","HALLMARK_APOPTOSIS","HALLMARK_HYPOXIA",
                  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_PROTEIN_SECRETION",
                  "HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_MYC_TARGETS_V1")
mods_k <- lapply(hm[keep_mods_pc], function(g) intersect(g, rownames(kdelta)))
mods_k$NEURONAL_CURATED <- intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                                       "GAP43","TUBB3","MAP2"), rownames(kdelta))
mods_k$PROLIF_CURATED <- mods$PROLIF
mods_k <- mods_k[sapply(mods_k, length) >= 8]
Md <- t(sapply(mods_k, function(g) colMeans(kdelta[g, , drop = FALSE], na.rm = TRUE)))
Md <- Md[, colSums(is.na(Md)) == 0]
pc <- prcomp(t(Md), scale. = TRUE)
pc2_load <- pc$rotation[, 2]

# GeoMx module scores per ROI (score the PC2 module set, keyed by full names)
score_roi <- t(sapply(names(pc2_load), function(mn) {
  g <- intersect(mods_k[[mn]], rownames(lnorm))
  colMeans(lnorm[g, , drop = FALSE], na.rm = TRUE)
}))
rownames(score_roi) <- names(pc2_load)
colnames(score_roi) <- colnames(lnorm)
pc2_roi <- as.numeric(t(pc2_load) %*% score_roi)
names(pc2_roi) <- colnames(score_roi)
# also the frozen 7 modules for display
score7 <- t(sapply(mods, function(g) colMeans(lnorm[g, , drop = FALSE], na.rm = TRUE)))

prox_scores <- t(sapply(proxies, function(g) colMeans(lnorm[g, , drop = FALSE], na.rm = TRUE)))
rownames(prox_scores) <- names(proxies)

# ---------- patient-level analysis ----------
res <- roi_meta
res$PC2 <- pc2_roi[res$dcc_stub]
for (mn in names(mods)) res[[paste0("mod_", mn)]] <- score7[mn, res$dcc_stub]
for (pn in names(proxies)) res[[paste0("px_", pn)]] <- prox_scores[pn, res$dcc_stub]

sink("results/geomx/spatial_second_axis.txt", split = TRUE)
cat("=== GeoMx spatial localization of the second axis ===\n")
cat("ROIs:", nrow(res), "| patients:", length(unique(res$patient)), "\n")
cat("stains: CD64 (myeloid-selected) vs GFAP (neural-selected)\n\n")

# patient-level: mean ROI score per patient x timepoint x stain
pt <- res %>% group_by(patient, timepoint, stain) %>%
  summarise(PC2 = mean(PC2), IL2STAT5 = mean(mod_IL2STAT5), ROS = mean(mod_ROS),
            EMT = mean(mod_EMT), Microglia = mean(px_Microglia),
            Endothelial = mean(px_Endothelial), .groups = "drop")

# 1. stain contrast within patients (CD64 vs GFAP ROIs)
wide <- pt %>% tidyr::pivot_wider(names_from = stain, values_from = c(PC2, IL2STAT5, ROS, EMT, Microglia, Endothelial))
wide <- wide[!is.na(wide$PC2_CD64) & !is.na(wide$PC2_GFAP), ]
cat("patients with both CD64 and GFAP ROIs:", nrow(wide), "\n")
cat("\n1. CD64 (myeloid) vs GFAP (neural) ROI contrast (paired within patient):\n")
for (v in c("PC2", "IL2STAT5", "ROS", "EMT", "Microglia", "Endothelial")) {
  d <- wide[[paste0(v, "_CD64")]] - wide[[paste0(v, "_GFAP")]]
  wt <- suppressWarnings(wilcox.test(d, mu = 0))
  cat("  ", v, ": median CD64-GFAP =", round(median(d), 2), " p =", signif(wt$p.value, 3), "\n")
}

# 2. primary vs recurrent (per stain, patient-level)
cat("\n2. primary vs recurrent per stain (patient-level):\n")
for (st in c("CD64", "GFAP")) {
  pt_st <- pt %>% filter(stain == st)
  for (v in c("PC2", "IL2STAT5", "ROS", "EMT")) {
    vp <- pt_st[[v]][pt_st$timepoint == "Primary"]
    vr <- pt_st[[v]][pt_st$timepoint == "Recurrent"]
    pp <- pt_st %>% dplyr::select(patient, timepoint, all_of(v)) %>%
      tidyr::pivot_wider(names_from = timepoint, values_from = all_of(v), values_fn = mean)
    pp <- pp[!is.na(pp$Primary) & !is.na(pp$Recurrent), ]
    wt_u <- suppressWarnings(wilcox.test(vr, vp))
    line <- paste0("  ", st, " ", v, ": primary n =", length(vp), " med ", round(median(vp), 2),
                   " | recurrent n =", length(vr), " med ", round(median(vr), 2),
                   " | unpaired p =", signif(wt_u$p.value, 3))
    if (nrow(pp) >= 3)
      line <- paste0(line, " | paired n =", nrow(pp), " delta =", round(median(pp$Recurrent - pp$Primary), 2),
                     " p =", signif(wilcox.test(pp$Recurrent - pp$Primary, mu = 0)$p.value, 3))
    cat(line, "\n")
  }
}

# 3. does PC2 co-localize with microglia/endothelial signal?
cat("\n3. PC2 vs compartment scores across ROIs:\n")
for (pn in c("px_Microglia", "px_Endothelial")) {
  ct <- suppressWarnings(cor.test(res$PC2, res[[pn]], method = "spearman"))
  cat("  PC2 ~", pn, ": rho =", round(ct$estimate, 2), " p =", signif(ct$p.value, 3), "\n")
}
sink()

# figures
p1 <- ggplot(pt, aes(stain, PC2, fill = stain)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = patient), alpha = 0.3, color = "grey60") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("CD64" = "#7B4EA3", "GFAP" = "#3D9B8C")) +
  facet_wrap(~ timepoint) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(y = "Frozen KNCC PC2 axis score (ROI level)", x = NULL,
       title = "Second axis localizes to CD64-selected (myeloid-rich) ROIs")
ggsave("figures_kncc/geomx_pc2_by_stain.pdf", p1, width = 7, height = 5)

p2 <- ggplot(res, aes(px_Microglia, PC2, color = stain)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("CD64" = "#7B4EA3", "GFAP" = "#3D9B8C")) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30") +
  theme_bw(base_size = 12) +
  labs(x = "Microglia score (ROI)", y = "PC2 axis score (ROI)",
       title = "PC2 axis tracks microglial signal across ROIs")
ggsave("figures_kncc/geomx_pc2_microglia.pdf", p2, width = 6, height = 5)

saveRDS(list(lnorm = lnorm, roi_meta = res, pt = pt), "results/geomx/geomx_bundle.rds")
cat("Done.\n")
