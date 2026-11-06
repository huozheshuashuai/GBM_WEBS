# =============================================================================
# 17b. STRING PPI via REST API (avoids STRINGdb download issues)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(igraph); library(ggplot2); library(dplyr) })

de <- read.csv("results/cptac_limma_results.csv")
splice_genes <- de$protein[de$is_spliceosome & !is.na(de$adj.P.Val) & de$adj.P.Val < 0.05]
WRITERS <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP")
ERASERS <- c("HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")
genes <- unique(c(WRITERS, ERASERS, splice_genes))
cat("Network input:", length(genes), "proteins\n")

# STRING REST API: network by identifiers (batch POST)
fetch_string <- function(gene_vec) {
  url <- "https://string-db.org/api/tsv/network"
  body <- list(identifiers = paste(gene_vec, collapse = "%0d"),
               species = 9606, required_score = 400,
               network_type = "functional")
  resp <- httr::POST(url, body = body)
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  read.delim(text = txt, stringsAsFactors = FALSE)
}

cat("Querying STRING API...\n")
net <- fetch_string(genes)
write.csv(net, "results/ppi_string_edges.csv", row.names = FALSE)
cat("Edges:", nrow(net), "\n")

# build graph
edges <- unique(data.frame(from = net$preferredName_A, to = net$preferredName_B,
                           score = net$score))
g <- graph_from_data_frame(edges, directed = FALSE)
V(g)$class <- ifelse(V(g)$name %in% WRITERS, "Writer",
              ifelse(V(g)$name %in% ERASERS, "Eraser", "Spliceosome"))

# topology
topo <- data.frame(
  gene = V(g)$name, class = V(g)$class,
  degree = degree(g),
  betweenness = betweenness(g),
  clustering = transitivity(g, type = "local")
) %>% arrange(desc(betweenness))
write.csv(topo, "results/ppi_topology.csv", row.names = FALSE)

# enrichment: observed edges vs expected (STRING ppi enrichment via API)
enrich_url <- "https://string-db.org/api/tsv/ppi_enrichment"
resp <- httr::POST(enrich_url, body = list(identifiers = paste(genes, collapse = "%0d"),
                                           species = 9606))
enr <- httr::content(resp, as = "text", encoding = "UTF-8")
cat("PPI enrichment:\n"); cat(enr, "\n")
writeLines(enr, "results/ppi_enrichment.txt")

# summary
sink("results/ppi_summary.txt", split = TRUE)
cat("=== STRING PPI network ===\n")
cat("Nodes:", vcount(g), " Edges:", ecount(g), "\n")
cat("Density:", round(edge_density(g), 4), "\n")
cat("Global clustering:", round(transitivity(g), 3), "\n")
cat("Components:", components(g)$no, "\n\n")
cat("Top 10 hubs by betweenness:\n"); print(head(topo, 10))
sink()

# plot
set.seed(42)
cols <- c("Writer" = "#d62728", "Eraser" = "#1f77b4", "Spliceosome" = "#9467bd")
pdf("figures_final/fig_ppi_network.pdf", width = 10, height = 9)
plot(g, vertex.color = cols[V(g)$class],
     vertex.size = scales::rescale(degree(g), to = c(4, 14)),
     vertex.label = V(g)$name, vertex.label.cex = 0.55, vertex.label.color = "black",
     edge.color = "grey75", main = sprintf("Writer-Eraser-Spliceosome PPI (%d nodes, %d edges)",
                                           vcount(g), ecount(g)))
legend("bottomleft", legend = names(cols), col = cols, pch = 19, bty = "n")
dev.off()
cat("Done.\n")
