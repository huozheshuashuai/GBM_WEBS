# =============================================================================
# 34. Fig 10: external validation (GSE131928 WEBS by Neftel state + platform)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(UCell); library(ggplot2); library(dplyr); library(patchwork) })
set.seed(42)
dir.create("figures_v2", showWarnings = FALSE)

WRITERS <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP")
ERASERS <- c("HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")
NEFTEL <- list(
  MES = c("CD44","VIM","S100A11","ANXA1","LGALS1","BCL2A1","CD83","CXCL8","AREG","CD74"),
  AC  = c("APOE","GFAP","S100B","HEPN1","GJA1","SPARCL1","CLU","AGT","FABP7","SLC1A3"),
  NPC = c("SOX4","STMN1","STMN2","BEX1","SOX11","DLL3","MARCKSL1","GAP43","INA","NNAT"),
  OPC = c("SOX10","PLP1","BCAS1","GPR17","PTPRZ1","CA2","MYRF","MAG","MOG","APOD")
)

read_tpm <- function(f) {
  m <- read.delim(gzfile(f), check.names = FALSE)
  g <- m[[1]]; m <- as.matrix(m[, -1]); rownames(m) <- make.unique(as.character(g))
  log2(m + 1)
}
ucell_diff <- function(mat) {
  w <- ScoreSignatures_UCell(mat, features = list(W = intersect(WRITERS, rownames(mat))))
  e <- ScoreSignatures_UCell(mat, features = list(E = intersect(ERASERS, rownames(mat))))
  wv <- if (nrow(w) == ncol(mat)) w[, 1] else w[1, ]
  ev <- if (nrow(e) == ncol(mat)) e[, 1] else e[1, ]
  wv - ev
}

dfs <- list()
for (f in c("GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz",
            "GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv.gz")) {
  platform <- ifelse(grepl("Smartseq2", f), "Smart-seq2", "10x")
  mat <- read_tpm(file.path("data/raw/GSE131928/extracted", f))
  webs <- ucell_diff(mat)
  mod_scores <- sapply(NEFTEL, function(gs) colMeans(mat[intersect(gs, rownames(mat)), , drop = FALSE]))
  state <- colnames(mod_scores)[max.col(mod_scores)]
  pts <- sub("^MGH", "", sub("[-_].*$", "", colnames(mat)))
  dfs[[platform]] <- data.frame(WEBS = webs, state = state, patient = pts, platform = platform)
}

all_df <- do.call(rbind, dfs)
all_df$state <- factor(all_df$state, levels = c("MES","AC","NPC","OPC"))

pA <- ggplot(all_df, aes(x = state, y = WEBS, fill = state)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_wrap(~platform) +
  scale_fill_manual(values = c(MES="#8c564b", AC="#2ca02c", NPC="#9467bd", OPC="#17becf")) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(title = "A. WEBS by Neftel state (GSE131928)", x = NULL, y = "WEBS")

pw <- all_df %>% group_by(patient, platform) %>% summarise(WEBS = median(WEBS), .groups = "drop")
wide <- tidyr::pivot_wider(pw, names_from = platform, values_from = WEBS) %>%
  filter(!is.na(`Smart-seq2`), !is.na(`10x`))
ct <- cor.test(wide$`Smart-seq2`, wide$`10x`, method = "spearman")
pB <- ggplot(wide, aes(x = `Smart-seq2`, y = `10x`)) +
  geom_point(size = 3) + geom_smooth(method = "lm", se = TRUE, color = "red") +
  theme_bw(base_size = 12) +
  labs(title = sprintf("B. Cross-platform consistency (rho = %.2f, n = %d)",
                       ct$estimate, nrow(wide)),
       x = "Patient WEBS (Smart-seq2)", y = "Patient WEBS (10x)")

pdf("figures_v2/Fig10_external_validation.pdf", width = 12, height = 5.5)
print(pA | pB)
dev.off()
message("Fig10 done")
