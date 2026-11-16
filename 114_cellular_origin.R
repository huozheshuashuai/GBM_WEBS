# =============================================================================
# 114. Task 1 (V6): cellular origin of the second axis (ROS/IL2-STAT5/EMT)
#      GSE138794 scRNA-seq: module expression by cell type
#      (author cell-type annotations; logCPM; per-sample then across-sample means)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

# author cell-type annotation (space-separated)
ct <- read.table("E:/KIMI/data/GSE138794_scRNA_Seq_cell_types.txt.gz",
                 col.names = c("barcode", "celltype"), stringsAsFactors = FALSE)
cat("annotated cells:", nrow(ct), "\n")
print(table(ct$celltype))

# module definitions (frozen 7)
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f2 <- strsplit(l, "\t")[[1]]; f2[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
mods <- list(
  PROLIF  = c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
              "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"),
  NEURO   = c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2","GAP43","TUBB3","MAP2"),
  ROS     = hm$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,
  IL2STAT5 = hm$HALLMARK_IL2_STAT5_SIGNALING,
  HYPOXIA = hm$HALLMARK_HYPOXIA,
  EMT     = hm$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  OXPHOS  = hm$HALLMARK_OXIDATIVE_PHOSPHORYLATION
)

# per-sample 10x processing
tenx_dirs <- list.dirs("data/raw/GSE138794/tenx", full.names = TRUE, recursive = FALSE)
cat("samples:", length(tenx_dirs), "\n")

library(Seurat)
ct_by_sample <- split(ct$celltype, sub("_.*", "", ct$barcode))

score_sample <- function(dir) {
  sid <- sub("GSM[0-9]*_", "", basename(dir))
  annot <- ct_by_sample[[sid]]
  if (is.null(annot)) return(NULL)
  counts <- tryCatch(Read10X(dir), error = function(e) NULL)
  if (is.null(counts)) return(NULL)
  bc <- paste0(sid, "_", sub("-[0-9]+$", "", colnames(counts)))
  keep <- bc %in% ct$barcode
  counts <- counts[, keep, drop = FALSE]
  celltypes <- ct$celltype[match(bc[keep], ct$barcode)]
  # logCPM
  lib <- colSums(counts)
  lcpm <- t(log2(t(counts) / (lib + 1) * 1e6 + 0.5))
  out <- lapply(names(mods), function(mn) {
    g <- intersect(mods[[mn]], rownames(lcpm))
    if (length(g) < 5) return(NULL)
    tapply(colMeans(lcpm[g, , drop = FALSE]), celltypes, mean)
  })
  names(out) <- names(mods)
  out
}

all_scores <- lapply(tenx_dirs, score_sample)
all_scores <- all_scores[!sapply(all_scores, is.null)]
cat("samples scored:", length(all_scores), "\n")

# aggregate: module x celltype mean across samples
ctypes <- sort(unique(unlist(lapply(all_scores, function(s) unlist(lapply(s, names))))))
agg <- lapply(names(mods), function(mn) {
  rows <- lapply(all_scores, function(s) s[[mn]])
  rows <- rows[!sapply(rows, is.null)]
  mat <- sapply(rows, function(v) v[ctypes])
  rownames(mat) <- ctypes
  rowMeans(mat, na.rm = TRUE)
})
names(agg) <- names(mods)
agg_df <- do.call(rbind, lapply(names(agg), function(mn)
  data.frame(module = mn, celltype = ctypes, score = agg[[mn]])))

sink("results/kncc/cellular_origin.txt", split = TRUE)
cat("=== cellular origin of second-axis modules (GSE138794 scRNA) ===\n\n")
for (mn in names(mods)) {
  d <- agg_df[agg_df$module == mn, ]
  d <- d[order(-d$score), ]
  cat(mn, ": top =", as.character(d$celltype[1]), round(d$score[1], 2),
      "| second =", as.character(d$celltype[2]), round(d$score[2], 2),
      "| ratio top/neoplastic =",
      round(d$score[d$celltype != "Neoplastic_cell"][1] / d$score[d$celltype == "Neoplastic_cell"], 2), "\n")
}
sink()

agg_df$module <- factor(agg_df$module, levels = names(mods))
p <- ggplot(agg_df, aes(celltype, module, fill = score)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#7B4EA3") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, fill = "Mean logCPM\nmodule score",
       title = "Cellular origin of recurrence-axis modules (GSE138794)")
ggsave("figures_kncc/cellular_origin.pdf", p, width = 7, height = 5)
write.csv(agg_df, "results/kncc/cellular_origin.csv", row.names = FALSE)
cat("Done.\n")
