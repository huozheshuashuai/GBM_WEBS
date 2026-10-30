# =============================================================================
# 22. CellChat v2: LLS vs LHS communication comparison
# =============================================================================
# Paper anchor: MIF -> CD74/CXCR4 macrophage-to-malignant signaling
# New: systematic LLS vs LHS differential communication (incl. LGALS9-HAVCR2,
#      SPP1-CD44), druggable pathway flagging
# Input: data/processed/GSE138794_final.rds
# Output: results/cellchat_*.csv, figures/fig_cellchat_*.pdf
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(Seurat); library(CellChat); library(dplyr); library(ggplot2) })
set.seed(42)

obj <- readRDS("data/processed/GSE138794_final.rds")

# unified cell-type label: prefer original annotation where available
obj$cc_label <- obj$cell_type
has_orig <- !is.na(obj$original_type)
obj$cc_label[has_orig] <- as.character(obj$original_type[has_orig])
obj$cc_label[obj$cc_label == "Neoplastic_cell"] <- "Malignant"
obj$cc_label[obj$cc_label == "Myeloid"] <- "Myeloid"
obj$cc_label[obj$cc_label == "Macrophage_Microglia"] <- "Myeloid"
obj$cc_label[obj$cc_label == "Endothelial_cell"] <- "Endothelial"
obj$cc_label[is.na(obj$cc_label)] <- "Unknown"

run_cellchat <- function(sub_obj, label) {
  input <- GetAssayData(sub_obj, layer = "data")
  meta <- data.frame(labels = sub_obj$cc_label, row.names = colnames(sub_obj))
  cc <- createCellChat(object = input, meta = meta, group.by = "labels")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

lls <- subset(obj, subset = WEBS_group == "LLS")
lhs <- subset(obj, subset = WEBS_group == "LHS")
message("LLS cells: ", ncol(lls), " | LHS cells: ", ncol(lhs))

cc_lls <- run_cellchat(lls, "LLS")
cc_lhs <- run_cellchat(lhs, "LHS")

# Pathway-level information flow comparison
info_lls <- cc_lls@netP$prob %>% reshape2::melt() %>%
  group_by(Var3) %>% summarise(flow_LLS = sum(value, na.rm = TRUE))
info_lhs <- cc_lhs@netP$prob %>% reshape2::melt() %>%
  group_by(Var3) %>% summarise(flow_LHS = sum(value, na.rm = TRUE))
colnames(info_lls)[1] <- colnames(info_lhs)[1] <- "pathway"
cmp <- full_join(info_lls, info_lhs, by = "pathway") %>%
  mutate(across(starts_with("flow"), ~ ifelse(is.na(.), 0, .)),
         delta = flow_LHS - flow_LLS,
         direction = ifelse(delta > 0, "Up in LHS", "Up in LLS")) %>%
  arrange(desc(delta))
write.csv(cmp, "results/cellchat_pathway_comparison.csv", row.names = FALSE)

# MIF pathway detail
df_mif <- tryCatch(subsetCommunication(cc_lhs, signaling = "MIF"),
                   error = function(e) { message("MIF not significant in LHS: ", e$message); NULL })
if (!is.null(df_mif)) write.csv(df_mif, "results/cellchat_MIF_LHS.csv", row.names = FALSE)

# Druggable flag
druggable <- c("MIF", "SPP1", "CD74", "CXCR4", "LGALS9", "HAVCR2", "VEGF", "PDGF")
cmp$druggable <- cmp$pathway %in% druggable

sink("results/cellchat_summary.txt")
cat("=== CellChat LLS vs LHS ===\n\n")
cat("LLS cells:", ncol(lls), "| LHS cells:", ncol(lhs), "\n\n")
cat("Top 10 pathways up in LHS:\n"); print(head(cmp, 10))
cat("\nTop 10 pathways up in LLS:\n"); print(tail(cmp, 10))
cat("\nDruggable pathways up in LHS:\n")
print(cmp[cmp$druggable & cmp$delta > 0, ])
sink()
cat(readLines("results/cellchat_summary.txt"), sep = "\n")

# merge & compare plot
merged <- tryCatch(mergeCellChat(list(LLS = cc_lls, LHS = cc_lhs), add.names = TRUE),
                   error = function(e) NULL)
dir.create("figures", showWarnings = FALSE)
pdf("figures/fig_cellchat_compare.pdf", width = 10, height = 6)
if (!is.null(merged)) {
  try(print(compareInteractions(merged, show.legend = FALSE, group = c(1, 2))),
      silent = TRUE)
  try(print(rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)),
      silent = TRUE)
}
dev.off()

saveRDS(list(LLS = cc_lls, LHS = cc_lhs), "results/cellchat_objects.rds")
message("Done.")
