# =============================================================================
# 17. PPI Network: Writers - Erasers - Spliceosome (STRING)
# =============================================================================
# Builds a STRING interaction network from writer/eraser enzymes plus the
# differentially abundant spliceosome proteins, and reports topology metrics.
# Input:  results/cptac_limma_results.csv (protein, logFC, adj.P.Val)
#         produced by 08_cptac_proteomics.R driver
# Output: results/ppi_topology.csv, figures/figS_ppi.pdf
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(STRINGdb)
  library(igraph)
  library(ggplot2)
  library(dplyr)
})

WRITER_GENES <- c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP")
ERASER_GENES <- c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")

# ---- Spliceosome proteins: significant from CPTAC limma (fallback: core list) ----
limma_file <- "results/cptac_limma_results.csv"
if (file.exists(limma_file)) {
  lim <- read.csv(limma_file)
  splice_genes <- lim$protein[lim$is_spliceosome & lim$adj.P.Val < 0.05]
} else {
  message("CPTAC limma results not found; using core spliceosome list from 08_cptac_proteomics.R")
  source("08_cptac_proteomics.R")
  splice_genes <- get_spliceosome_genes()
}

genes <- unique(c(WRITER_GENES, ERASER_GENES, splice_genes))
message("Network input: ", length(genes), " proteins")

# ---- STRING ----
string_db <- STRINGdb$new(version = "12.0", species = 9606,
                          score_threshold = 400, network_type = "full",
                          input_directory = "data/raw/STRING")
mapped <- string_db$map(data.frame(gene = genes), "gene", removeUnmappedRows = TRUE)

# interaction table
interactions <- string_db$get_interactions(mapped$STRING_id)
g <- graph_from_data_frame(interactions, directed = FALSE, vertices = mapped$STRING_id)
V(g)$gene <- mapped$gene[match(V(g)$name, mapped$STRING_id)]
V(g)$class <- ifelse(V(g)$gene %in% WRITER_GENES, "Writer",
              ifelse(V(g)$gene %in% ERASER_GENES, "Eraser", "Spliceosome"))

# ---- Topology metrics ----
topo <- data.frame(
  gene = V(g)$gene,
  class = V(g)$class,
  degree = degree(g),
  betweenness = betweenness(g),
  clustering = transitivity(g, type = "local")
) %>% arrange(desc(betweenness))

dir.create("results", showWarnings = FALSE)
write.csv(topo, "results/ppi_topology.csv", row.names = FALSE)

sink("results/ppi_summary.txt")
cat("=== PPI network summary (STRING v12, score >= 400) ===\n")
cat("Nodes:", vcount(g), " Edges:", ecount(g), "\n")
cat("Density:", round(edge_density(g), 4), "\n")
cat("Global clustering coefficient:", round(transitivity(g), 4), "\n")
cat("Components:", components(g)$no, "\n\n")
cat("Observed edges vs random expectation (STRING enrichment):\n")
print(string_db$get_ppi_enrichment(mapped$STRING_id))
cat("\nTop 10 hubs by betweenness:\n")
print(head(topo, 10))
sink()

# ---- Plot ----
set.seed(42)
dir.create("figures", showWarnings = FALSE)
cols <- c("Writer" = "#d62728", "Eraser" = "#1f77b4", "Spliceosome" = "#ff7f0e")
pdf("figures/figS_ppi.pdf", width = 10, height = 9)
plot(g,
     vertex.color = cols[V(g)$class],
     vertex.size = scales::rescale(degree(g), to = c(4, 14)),
     vertex.label = V(g)$gene,
     vertex.label.cex = 0.55,
     vertex.label.color = "black",
     edge.color = "grey70",
     main = sprintf("Writer-Eraser-Spliceosome PPI network (%d nodes, %d edges)",
                    vcount(g), ecount(g)))
legend("bottomleft", legend = names(cols), col = cols, pch = 19, bty = "n")
dev.off()

message("Done. See results/ppi_topology.csv and figures/figS_ppi.pdf")
