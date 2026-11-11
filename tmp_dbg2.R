suppressMessages({ library(dplyr) })
source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
prot <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct", check.names = FALSE)
pm <- as.matrix(sapply(prot[, -1], as.numeric)); rownames(pm) <- prot[[1]]
rna <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct", check.names = FALSE)
rm_ <- as.matrix(sapply(rna[, -1], as.numeric)); rownames(rm_) <- rna[[1]]
genes <- intersect(rownames(pm), rownames(rm_))
webs <- read.csv("results/cptac_webs.csv")
samples <- intersect(colnames(pm), intersect(colnames(rm_), webs$sample))
cat("samples:", length(samples), "\n")
pm2 <- pm[genes, samples]; rm2 <- rm_[genes, samples]
cat("pm2 dim:", dim(pm2), " rm2 dim:", dim(rm2), "\n")
resid_m <- do.call(rbind, lapply(genes, function(g) resid(lm(pm2[g, ] ~ rm2[g, ]))))
cat("resid_m dim:", dim(resid_m), " class:", class(resid_m), "\n")
sp <- intersect(splice, genes)
cat("sp:", length(sp), "\n")
sr <- colMeans(resid_m[sp, , drop = FALSE], na.rm = TRUE)
cat("resid score len:", length(sr), "\n")
