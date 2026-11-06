# =============================================================================
# 57. PTM abundance correction: are acetylation/phosphorylation changes
#     independent of parent-protein abundance?
#     Corrected PTM = site intensity - parent protein abundance (log scale)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(limma); library(dplyr) })

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
webs <- read.csv("results/cptac_webs.csv")

# parent protein abundance
prot <- read.delim("E:/KIMI/data/HS_CPTAC_GBM_proteome_mssm_per_gene.cct", check.names = FALSE)
pmat <- as.matrix(sapply(prot[, -1], as.numeric)); rownames(pmat) <- prot[[1]]

run_correction <- function(site_file, label) {
  ptm <- read.delim(site_file, check.names = FALSE)
  sites <- ptm[[1]]
  smat <- as.matrix(sapply(ptm[, -1], as.numeric)); rownames(smat) <- sites
  site_gene <- sub("-.*$", "", sites)

  common <- Reduce(intersect, list(colnames(smat), colnames(pmat), webs$sample))
  smat <- smat[, common]
  w <- webs[match(common, webs$sample), ]
  grp <- factor(w$WEBS_group, levels = c("LLS", "LDTS", "LHS"))

  # corrected = site - matched protein abundance
  ok <- !is.na(site_gene) & site_gene %in% rownames(pmat)
  smat_c <- smat[ok, ]
  site_gene_c <- site_gene[ok]
  prot_matched <- pmat[site_gene_c, common]
  smat_corr <- smat_c - prot_matched   # log-ratio to parent protein

  # limma LHS vs LLS on RAW and CORRECTED
  sel <- grp %in% c("LLS", "LHS")
  g_sel <- droplevels(grp[sel])
  des <- model.matrix(~ g_sel)

  run_limma <- function(mat) {
    keep <- rowMeans(!is.na(mat)) >= 0.5
    fit <- eBayes(lmFit(mat[keep, sel], des))
    tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    tt$site <- rownames(tt)
    tt$gene <- sub("-.*$", "", tt$site)
    tt$is_spliceosome <- tt$gene %in% splice
    tt
  }
  tt_raw <- local({
    keep <- rowMeans(!is.na(smat)) >= 0.5
    fit <- eBayes(lmFit(smat[keep, sel], des))
    tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    tt$site <- rownames(tt); tt$gene <- sub("-.*$", "", tt$site)
    tt$is_spliceosome <- tt$gene %in% splice; tt
  })
  tt_cor <- local({
    keep <- rowMeans(!is.na(smat_corr)) >= 0.5
    fit <- eBayes(lmFit(smat_corr[keep, sel], des))
    tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    tt$site <- rownames(tt); tt$gene <- sub("-.*$", "", tt$site)
    tt$is_spliceosome <- tt$gene %in% splice; tt
  })

  fisher_or <- function(tt) {
    sig <- tt$adj.P.Val < 0.05
    tab <- matrix(c(sum(sig & tt$is_spliceosome), sum(sig & !tt$is_spliceosome),
                    sum(!sig & tt$is_spliceosome), sum(!sig & !tt$is_spliceosome)), 2, byrow = TRUE)
    ft <- fisher.test(tab)
    c(OR = unname(ft$estimate), p = ft$p.value,
      n_sig = sum(sig), n_sig_splice = sum(sig & tt$is_spliceosome),
      n_splice = sum(tt$is_spliceosome))
  }
  r0 <- fisher_or(tt_raw); r1 <- fisher_or(tt_cor)

  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("RAW:       sig=%d, spliceosome sig=%d/%d, OR=%.2f p=%.3g\n",
              r0["n_sig"], r0["n_sig_splice"], r0["n_splice"], r0["OR"], r0["p"]))
  cat(sprintf("CORRECTED: sig=%d, spliceosome sig=%d/%d, OR=%.2f p=%.3g\n",
              r1["n_sig"], r1["n_sig_splice"], r1["n_splice"], r1["OR"], r1["p"]))

  # per-site coupling: acetylation ~ abundance correlation (key claim)
  if (label == "acetylome") {
    # does site acetylation still correlate with WEBS after correction?
    sp_sites <- tt_cor$site[tt_cor$is_spliceosome]
    sp_sites <- intersect(sp_sites, rownames(smat_corr))
    cors <- sapply(sp_sites, function(s) {
      x <- smat_corr[s, ]
      ok <- !is.na(x)
      if (sum(ok) < 30) return(NA)
      cor(x[ok], w$WEBS[ok], method = "spearman")
    })
    cat(sprintf("Spliceosome sites (corrected) correlated with WEBS: median rho=%.3f, FDR<0.05: %d/%d\n",
                median(cors, na.rm = TRUE),
                sum(p.adjust(sapply(sp_sites, function(s) {
                  x <- smat_corr[s, ]; ok <- !is.na(x)
                  if (sum(ok) < 30) return(NA)
                  suppressWarnings(cor.test(x[ok], w$WEBS[ok], method = "spearman"))$p.value
                }), method = "BH") < 0.05, na.rm = TRUE),
                length(sp_sites)))
  }
  write.csv(tt_cor, paste0("results/", label, "_abundance_corrected.csv"), row.names = FALSE)
  invisible(list(raw = r0, cor = r1))
}

sink("results/ptm_abundance_correction.txt", split = TRUE)
cat("=== PTM abundance correction (site - parent protein) ===\n")
a <- run_correction("E:/KIMI/data/HS_CPTAC_GBM_acetylome_mssm_site.cct", "acetylome")
p <- run_correction("E:/KIMI/data/HS_CPTAC_GBM_phosphoproteome_mssm_site.cct", "phosphoproteome")
sink()
cat("Done.\n")
