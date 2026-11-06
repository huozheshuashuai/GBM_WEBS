suppressMessages(library(SNFtool))
load_view <- function(f, collapse = FALSE) {
  m <- as.matrix(read.delim(f, row.names = 1, check.names = FALSE))
  m <- m[rowSums(!is.na(m)) > ncol(m) * 0.5, ]
  if (collapse) {
    g <- sub("[-_\\.][A-Za-z]*[0-9]+.*$", "", rownames(m))
    m <- rowsum(m, g, reorder = FALSE, na.rm = TRUE) /
         pmax(1, as.vector(table(g)[rownames(rowsum(m, g, reorder = FALSE))]))
  }
  m
}
topv <- function(m, n = 2000) m[order(apply(m, 1, var, na.rm = TRUE), decreasing = TRUE)[1:min(n, nrow(m))], ]
prot <- load_view("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct")
acet <- load_view("data/raw/CPTAC/HS_CPTAC_GBM_acetylome_mssm_site.cct", collapse = TRUE)
phos <- load_view("data/raw/CPTAC/HS_CPTAC_GBM_phosphoproteome_mssm_site.cct", collapse = TRUE)
rna <- load_view("data/raw/CPTAC/HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct")
views <- list(proteome = topv(prot), acetylome = topv(acet), phospho = topv(phos), rnaseq = topv(rna))
samples <- Reduce(intersect, lapply(views, colnames))
views <- lapply(views, function(m) { m <- m[, samples]; m[complete.cases(m), ] })
cat("samples:", length(samples), "\n")
cat("views dims:", sapply(views, function(x) paste(dim(x), collapse = "x")), "\n")
Ws <- lapply(views, function(m) { d <- dist(t(scale(t(m)))); affinityMatrix(as.matrix(d), K = 20, sigma = 0.5) })
cat("W dims:", sapply(Ws, function(x) paste(dim(x), collapse = "x")), "\n")
cat("W has NA:", sapply(Ws, function(x) sum(is.na(x))), "\n")
