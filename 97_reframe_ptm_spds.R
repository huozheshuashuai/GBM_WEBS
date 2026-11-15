# =============================================================================
# 97. WEBS-free PTM analysis: spliceosome acetylation/phosphorylation vs SPDS
#     (spliceosome protein-mRNA discordance score), raw + parent-protein-corrected
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

ax <- readRDS("results/reframe/cptac_discordance_axis.rds")
spds <- ax$spds   # named by sample

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()

proteome <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct",
                                 row.names = 1, check.names = FALSE))

analyze_ptm <- function(file, site_regex, label, min_n = 40) {
  cat("\n#####", label, "#####\n")
  M <- as.matrix(read.delim(file, row.names = 1, check.names = FALSE))
  site_ids <- rownames(M)
  site_gene <- sub(site_regex, "", site_ids)
  samples <- intersect(colnames(M), names(spds))
  cat("samples matched:", length(samples), "\n")
  M <- M[, samples]
  x <- spds[samples]
  ok <- rowSums(!is.na(M)) >= min_n
  M <- M[ok, ]; sg <- site_gene[ok]
  cat("sites tested:", nrow(M), "| spliceosome sites:", sum(sg %in% splice), "\n")

  # raw correlation with SPDS
  rr <- t(apply(M, 1, function(v) {
    keep <- !is.na(v)
    if (sum(keep) < min_n) return(c(NA, NA))
    ct <- suppressWarnings(cor.test(v[keep], x[keep], method = "spearman"))
    c(ct$estimate, ct$p.value)
  }))
  res <- data.frame(site = rownames(M), gene = sg, rho = rr[, 1], p = rr[, 2])
  res$q <- p.adjust(res$p, "BH")
  res$is_spliceosome <- res$gene %in% splice

  enrich <- function(df, qcol, rho_col, tag) {
    sigpos <- !is.na(df[[qcol]]) & df[[qcol]] < 0.05 & df[[rho_col]] > 0
    tab <- matrix(c(sum(sigpos & df$is_spliceosome), sum(!sigpos & df$is_spliceosome),
                    sum(sigpos & !df$is_spliceosome), sum(!sigpos & !df$is_spliceosome)), 2)
    ft <- fisher.test(tab)
    cat(tag, ": spliceosome sites positively SPDS-associated:",
        sum(sigpos & df$is_spliceosome), "/", sum(df$is_spliceosome),
        "| OR =", round(ft$estimate, 2), " p =", signif(ft$p.value, 3), "\n")
    ft
  }
  enrich(res, "q", "rho", "RAW")

  # parent-protein correction (site - parent protein, log scale)
  psam <- intersect(samples, colnames(proteome))
  pshare <- proteome[, psam]
  parent <- pshare[match(res$gene, rownames(pshare)), ]
  Mc <- M[, psam] - parent
  xc <- spds[psam]
  okc <- rowSums(!is.na(Mc)) >= min_n
  Mc <- Mc[okc, ]
  res_c <- res[okc, ]
  rr2 <- t(sapply(seq_len(nrow(Mc)), function(i) {
    v <- Mc[i, ]; keep <- !is.na(v)
    if (sum(keep) < min_n) return(c(NA, NA))
    ct <- suppressWarnings(cor.test(v[keep], xc[keep], method = "spearman"))
    c(ct$estimate, ct$p.value)
  }))
  res_c$rho_adj <- rr2[, 1]; res_c$p_adj <- rr2[, 2]
  res_c$q_adj <- p.adjust(res_c$p_adj, "BH")
  enrich(res_c, "q_adj", "rho_adj", "CORRECTED")

  # median rho of spliceosome sites
  cat("spliceosome median rho: raw", round(median(res$rho[res$is_spliceosome], na.rm = TRUE), 3),
      "| corrected", round(median(res_c$rho_adj[res_c$is_spliceosome], na.rm = TRUE), 3), "\n")
  write.csv(res_c, paste0("results/reframe/", tolower(gsub(" ", "_", label)), "_spds.csv"),
            row.names = FALSE)
  invisible(res_c)
}

acet <- analyze_ptm("data/raw/CPTAC/HS_CPTAC_GBM_acetylome_mssm_site.cct",
                    "[-_\\.][Kk][0-9]+.*$", "Acetylome")
phos <- analyze_ptm("data/raw/CPTAC/HS_CPTAC_GBM_phosphoproteome_mssm_site.cct",
                    "[-_\\.][STYsty][0-9]+.*$", "Phosphoproteome")
cat("Done.\n")
