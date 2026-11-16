setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(Seurat))
ct <- read.table("E:/KIMI/data/GSE138794_scRNA_Seq_cell_types.txt.gz",
                 col.names = c("barcode", "celltype"), stringsAsFactors = FALSE)
ct_by_sample <- split(ct$celltype, sub("_.*", "", ct$barcode))
dirs <- list.dirs("data/raw/GSE138794/tenx", full.names = TRUE, recursive = FALSE)
for (dir in dirs) {
  sid <- sub("GSM[0-9]*_", "", basename(dir))
  annot <- ct_by_sample[[sid]]
  if (is.null(annot)) { cat(sid, ": no annot\n"); next }
  counts <- tryCatch(Read10X(dir), error = function(e) { cat(sid, "Read10X error:", conditionMessage(e), "\n"); NULL })
  if (is.null(counts)) next
  bc <- paste0(sid, "_", sub("-[0-9]+$", "", colnames(counts)))
  keep <- bc %in% ct$barcode
  cat(sid, ": cells", ncol(counts), "| matched", sum(keep), "\n")
}
