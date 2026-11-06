# =============================================================================
# 49. CPTAC phosphoproteome: spliceosome phosphorylation & PTM crosstalk
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(limma); library(dplyr); library(ggplot2) })

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()

# --- load phosphosite matrix ---
ph <- read.delim("E:/KIMI/data/HS_CPTAC_GBM_phosphoproteome_mssm_site.cct",
                 check.names = FALSE)
sites <- ph[[1]]
pmat <- as.matrix(sapply(ph[, -1], as.numeric))
rownames(pmat) <- sites
site_gene <- sub("-.*$", "", sites)

# --- WEBS groups ---
webs <- read.csv("results/cptac_webs.csv")
common <- intersect(colnames(pmat), webs$sample)
pmat <- pmat[, common]
webs <- webs[match(common, webs$sample), ]
cat("Matched samples:", length(common), "\n")
grp <- factor(webs$WEBS_group, levels = c("LLS", "LDTS", "LHS"))
cat("Group sizes:", table(grp), "\n")

# keep sites quantified in >=50% samples
keep <- rowMeans(!is.na(pmat)) >= 0.5
pmat2 <- pmat[keep, ]
cat("Sites after filtering:", nrow(pmat2), "\n")

# --- limma: LHS vs LLS ---
sel <- grp %in% c("LLS", "LHS")
g_sel <- droplevels(grp[sel])
design <- model.matrix(~ g_sel)
fit <- lmFit(pmat2[, sel], design)
fit <- eBayes(fit)
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
tt$site <- rownames(tt)
tt$gene <- sub("-.*$", "", tt$site)
tt$is_spliceosome <- tt$gene %in% splice
write.csv(tt, "results/phospho_limma_LHS_vs_LLS.csv", row.names = FALSE)

sink("results/phospho_summary.txt", split = TRUE)
cat("=== CPTAC phosphoproteome: LHS vs LLS ===\n")
cat("Sites tested:", nrow(tt), "\n")
sig <- tt[tt$adj.P.Val < 0.05, ]
cat("Significant sites (FDR<0.05):", nrow(sig), "\n")
cat("Spliceosome sites tested:", sum(tt$is_spliceosome),
    " significant:", sum(sig$is_spliceosome), "\n")
# enrichment of spliceosome among significant (Fisher)
tab <- matrix(c(sum(sig$is_spliceosome), nrow(sig) - sum(sig$is_spliceosome),
                sum(tt$is_spliceosome) - sum(sig$is_spliceosome),
                nrow(tt) - nrow(sig) - sum(tt$is_spliceosome) + sum(sig$is_spliceosome)),
              nrow = 2, byrow = TRUE)
ft <- fisher.test(tab)
cat("Fisher spliceosome enrichment: OR =", round(ft$estimate, 2),
    " p =", signif(ft$p.value, 3), "\n\n")
cat("Top spliceosome phosphosites (by |t|):\n")
print(tt %>% filter(is_spliceosome) %>% arrange(desc(abs(t))) %>% head(15) %>%
        dplyr::select(site, logFC, t, adj.P.Val))

# --- SR protein phosphosites (splicing regulators) ---
sr_genes <- paste0("SRSF", 1:12)
sr_tt <- tt[tt$gene %in% sr_genes, ]
cat("\nSRSF family phosphosites:", nrow(sr_tt),
    " significant:", sum(sr_tt$adj.P.Val < 0.05), "\n")
print(sr_tt %>% arrange(adj.P.Val) %>% head(10) %>%
        dplyr::select(site, logFC, adj.P.Val))

# --- kinase abundance vs WEBS (SRPK/CLK from proteome) ---
prot <- read.delim("E:/KIMI/data/HS_CPTAC_GBM_proteome_mssm_per_gene.cct",
                   check.names = FALSE)
pg <- prot[[1]]
pm <- as.matrix(sapply(prot[, -1], as.numeric))
rownames(pm) <- pg
kinases <- c("SRPK1","SRPK2","CLK1","CLK2","CLK3","CLK4","CDK7","CDK9","CDK12","CDK13")
kinases <- intersect(kinases, rownames(pm))
kk <- intersect(colnames(pm), webs$sample)
cat("\nKinase abundance vs WEBS (Spearman):\n")
for (k in kinases) {
  ct <- suppressWarnings(cor.test(pm[k, kk], webs$WEBS[match(kk, webs$sample)],
                                  method = "spearman"))
  cat(sprintf("  %-7s rho = %+.3f  p = %.3g\n", k, ct$estimate, ct$p.value))
}

# --- PTM crosstalk: acetylation vs phosphorylation on same spliceosome proteins ---
ac <- read.delim("E:/KIMI/data/HS_CPTAC_GBM_acetylome_mssm_site.cct",
                 check.names = FALSE)
ac_gene <- sub("-.*$", "", ac[[1]])
splice_ac <- unique(ac_gene[ac_gene %in% splice])
splice_ph <- unique(tt$gene[tt$is_spliceosome])
both <- intersect(splice_ac, splice_ph)
cat("\nSpliceosome proteins with BOTH acetylation and phosphorylation data:",
    length(both), "\n")
# per-protein: mean |t| of acetyl sites vs phospho sites
ac_sites <- ac[[1]]
acmat <- as.matrix(sapply(ac[, -1], as.numeric))
rownames(acmat) <- ac_sites
ck <- intersect(colnames(acmat), webs$sample)
gx <- factor(webs$WEBS_group[match(ck, webs$sample)], levels = c("LLS","LDTS","LHS"))
s2 <- gx %in% c("LLS","LHS")
g_s2 <- droplevels(gx[s2])
design2 <- model.matrix(~ g_s2)
fit2 <- eBayes(lmFit(acmat[rowMeans(!is.na(acmat)) >= 0.5, ck[s2]], design2))
tt2 <- topTable(fit2, coef = 2, number = Inf, sort.by = "none")
tt2$gene <- sub("-.*$", "", rownames(tt2))
xt <- data.frame(gene = both,
                 acetyl_t = sapply(both, function(g) mean(abs(tt2$t[tt2$gene == g]), na.rm = TRUE)),
                 phospho_t = sapply(both, function(g) mean(abs(tt$t[tt$gene == g]), na.rm = TRUE)))
xt <- xt[is.finite(xt$acetyl_t) & is.finite(xt$phospho_t), ]
ct <- cor.test(xt$acetyl_t, xt$phospho_t, method = "spearman")
cat("Crosstalk: per-protein |t| acetylation vs phosphorylation: rho =",
    round(ct$estimate, 3), " p =", signif(ct$p.value, 3),
    " (n =", nrow(xt), "proteins)\n")
sink()
write.csv(xt, "results/ptm_crosstalk_spliceosome.csv", row.names = FALSE)

# --- figure: volcano colored by spliceosome ---
tt$sig_flag <- ifelse(tt$adj.P.Val < 0.05 & tt$is_spliceosome, "Spliceosome sig",
               ifelse(tt$adj.P.Val < 0.05, "Other sig",
               ifelse(tt$is_spliceosome, "Spliceosome NS", "NS")))
p <- ggplot(tt, aes(logFC, -log10(P.Value), color = sig_flag)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Spliceosome sig" = "#7B4EA3", "Other sig" = "#D4A017",
                                "Spliceosome NS" = "#3D9B8C", "NS" = "grey85")) +
  theme_bw(base_size = 12) +
  labs(title = "Phosphosite changes LHS vs LLS (CPTAC)",
       subtitle = paste0(sum(tt$is_spliceosome & tt$adj.P.Val < 0.05),
                         " significant spliceosome phosphosites"),
       x = "logFC (LHS vs LLS)", y = "-log10(p)", color = NULL)
ggsave("figures_final/fig_phospho_volcano.pdf", p, width = 7, height = 5)
cat("Done.\n")
