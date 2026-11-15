# =============================================================================
# 90. KNCC module 5: spliceosome phosphosites in paired T1->T2 (REVISED)
#     raw site deltas + TWO parent-protein correction methods:
#       A) delta subtraction (site_delta - parent_delta), as in CPTAC analysis
#       B) regression residual (site_delta ~ parent_delta across patients)
#     plus site~parent tracking diagnostics and spliceosome enrichment
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2); library(org.Hs.eg.db); library(AnnotationDbi) })
set.seed(42)

base <- "E:/KIMI/data/PDC000515/"
f <- paste0(base, "KNCC_Glioblastoma_Evolution_Phosphoproteome.phosphosite.tmt11.tsv")
hdr <- strsplit(readLines(f, n = 1), "\t")[[1]]
log_idx <- which(grepl("^Log KNCC_GBM", hdr))
gene_col <- which(hdr == "Gene")
cc <- rep("NULL", length(hdr)); cc[1] <- "character"; cc[gene_col] <- "character"; cc[log_idx] <- "numeric"
ph <- read.delim(f, sep = "\t", header = TRUE, colClasses = cc, check.names = FALSE)
cat("phosphosite rows:", nrow(ph), "\n")
site_gene <- ph$Gene

colnames(ph)[-1] <- sub("\\.[0-9]+$", "", sub(":.*", "", sub("^Log ", "", colnames(ph)[-1])))
num_cols <- which(grepl("^KNCC_GBM", colnames(ph)))
m <- as.matrix(ph[, num_cols]); cn <- colnames(m)
mu <- matrix(NA_real_, nrow(m), length(unique(cn)), dimnames = list(NULL, unique(cn)))
for (s in unique(cn)) { idx <- which(cn == s); mu[, s] <- if (length(idx) == 1) m[, idx] else rowMeans(m[, idx, drop = FALSE], na.rm = TRUE) }
phm <- mu; rownames(phm) <- ph[[1]]
cat("phospho matrix:", nrow(phm), "sites x", ncol(phm), "samples\n")

# pairing
samples <- colnames(phm)
patient <- sub("_(T[12])$", "", samples); tp <- sub(".*_", "", samples)
ptab <- table(patient, tp)
paired_pts <- rownames(ptab)[ptab[, "T1"] > 0 & ptab[, "T2"] > 0]
t1 <- samples[tp == "T1"][match(paired_pts, patient[tp == "T1"])]
t2 <- samples[tp == "T2"][match(paired_pts, patient[tp == "T2"])]
dph <- phm[, t2] - phm[, t1]
colnames(dph) <- paired_pts
cat("phospho patients:", length(unique(patient)), "| paired T1+T2:", length(paired_pts), "\n")

# gene symbols come with the file (no ENSP mapping needed)
gene <- site_gene
cat("sites with gene annotation:", sum(!is.na(gene) & gene != ""), "/", length(gene), "\n")

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
is_sp <- gene %in% splice
is_sp[is.na(is_sp)] <- FALSE
cat("spliceosome sites:", sum(is_sp), "\n")

# parent-protein deltas (proteome paired patients)
obj <- readRDS("results/kncc/proteome_analysis_ready.rds")
prot_delta <- obj$delta
common_pts <- intersect(colnames(prot_delta), colnames(dph))
cat("patients paired in BOTH proteome and phospho:", length(common_pts), "\n")
dph_c <- dph[, common_pts]
parent <- prot_delta[, common_pts]
parent_delta <- parent[match(gene, rownames(parent)), ]
has_parent <- !is.na(gene) & gene %in% rownames(parent)
cat("sites with parent protein quantified:", sum(has_parent), "/", nrow(dph_c), "\n")

site_test <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) < 25) return(c(med = NA_real_, p = NA_real_))
  c(med = median(v), p = suppressWarnings(wilcox.test(v, mu = 0))$p.value)
}

# ---- 1. raw site deltas ----
rr <- t(apply(dph, 1, site_test))
site_df <- data.frame(site = rownames(dph), gene = gene, is_spliceosome = is_sp,
                      med = rr[, 1], p = rr[, 2])
site_df$q <- p.adjust(site_df$p, "BH")

# ---- 2a. correction A: delta subtraction ----
dph_sub <- dph_c - parent_delta
rrA <- t(apply(dph_sub, 1, site_test))
site_df$med_sub <- NA; site_df$p_sub <- NA
site_df$med_sub[match(rownames(dph_sub), site_df$site)] <- rrA[, 1]
site_df$p_sub[match(rownames(dph_sub), site_df$site)] <- rrA[, 2]
site_df$q_sub <- p.adjust(site_df$p_sub, "BH")

# ---- 2b. correction B: regression residual across patients ----
resid_reg <- matrix(NA_real_, nrow(dph_c), ncol(dph_c), dimnames = dimnames(dph_c))
for (i in seq_len(nrow(dph_c))) {
  ok <- !is.na(dph_c[i, ]) & !is.na(parent_delta[i, ])
  if (sum(ok) < 25) next
  fit <- tryCatch(lm(dph_c[i, ok] ~ parent_delta[i, ok]), error = function(e) NULL)
  if (!is.null(fit)) resid_reg[i, ok] <- resid(fit)
}
rrB <- t(apply(resid_reg, 1, site_test))
site_df$med_reg <- NA; site_df$p_reg <- NA
site_df$med_reg[match(rownames(resid_reg), site_df$site)] <- rrB[, 1]
site_df$p_reg[match(rownames(resid_reg), site_df$site)] <- rrB[, 2]
site_df$q_reg <- p.adjust(site_df$p_reg, "BH")

# ---- diagnostics: how strongly do sites track their parent protein? ----
track <- sapply(which(has_parent)[1:min(3000, sum(has_parent))], function(i) {
  ok <- !is.na(dph_c[i, ]) & !is.na(parent_delta[i, ])
  suppressWarnings(cor(dph_c[i, ok], parent_delta[i, ok], method = "spearman"))
})
cat("site-parent tracking (Spearman, 3000 sampled sites): median rho =",
    round(median(track, na.rm = TRUE), 3), "\n")

# ---- spliceosome summary ----
spsum <- function(med, q, label) {
  s <- site_df$is_spliceosome & !is.na(med)
  sig <- s & q < 0.05
  cat(label, ": spliceosome n =", sum(s),
      " | median =", round(median(med[s]), 3),
      " | FDR<0.05:", sum(sig),
      " (down:", sum(sig & med < 0), " up:", sum(sig & med > 0), ")\n")
  wt <- wilcox.test(med[s], med[!s & !is.na(med)])
  cat("   spliceosome vs genome: p =", signif(wt$p.value, 3), "\n")
  # enrichment among significant sites
  ft <- fisher.test(matrix(c(sum(sig), sum(s & !(q < 0.05)),
                             sum(!site_df$is_spliceosome & !is.na(med) & q < 0.05, na.rm = TRUE),
                             sum(!site_df$is_spliceosome & !is.na(med) & !(q < 0.05), na.rm = TRUE)), 2, byrow = TRUE))
  cat("   enrichment OR =", round(ft$estimate, 2), " p =", signif(ft$p.value, 3), "\n")
}
sink("results/kncc/module5_phospho.txt", split = TRUE)
cat("=== Module 5: paired phosphosite analysis (KNCC) ===\n")
cat("sites:", nrow(site_df), "| spliceosome sites:", sum(is_sp),
    "| paired phospho patients:", length(paired_pts),
    "| with paired proteome:", length(common_pts), "\n")
cat("site-parent tracking median rho:", round(median(track, na.rm = TRUE), 3), "\n\n")
spsum(site_df$med, site_df$q, "RAW          ")
spsum(site_df$med_sub, site_df$q_sub, "SUBTRACTION  ")
spsum(site_df$med_reg, site_df$q_reg, "REGRESSION   ")

# SRSF family
srsf <- grepl("^SRSF", site_df$gene)
srsf[is.na(srsf)] <- FALSE
cat("\nSRSF sites n =", sum(srsf),
    " | raw median", round(median(site_df$med[srsf], na.rm = TRUE), 3),
    " | reg-corrected median", round(median(site_df$med_reg[srsf], na.rm = TRUE), 3),
    " | reg FDR<0.05:", sum(site_df$q_reg[srsf] < 0.05, na.rm = TRUE),
    " (down:", sum(site_df$med_reg[srsf] < 0 & site_df$q_reg[srsf] < 0.05, na.rm = TRUE), ")\n")

# top corrected spliceosome sites
top <- site_df %>% filter(is_spliceosome, !is.na(q_reg)) %>% arrange(p_reg) %>% head(15)
cat("\ntop 15 spliceosome sites by corrected p:\n")
print(top[, c("site", "gene", "med", "med_reg", "q_reg")], row.names = FALSE)
sink()

write.csv(site_df, "results/kncc/phosphosite_paired.csv", row.names = FALSE)

# figure: raw vs corrected (regression) for spliceosome sites
plot_df <- site_df %>% filter(!is.na(med), !is.na(med_reg))
p <- ggplot(plot_df, aes(med, med_reg)) +
  geom_point(color = "grey85", size = 0.4) +
  geom_point(data = plot_df %>% filter(is_spliceosome), color = "#7B4EA3", size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey60") + geom_vline(xintercept = 0, color = "grey60") +
  theme_bw(base_size = 11) +
  labs(x = "Raw site delta (T2 - T1)", y = "Abundance-corrected site delta (regression)",
       title = "Phosphosite changes: raw vs parent-protein-corrected",
       subtitle = sprintf("spliceosome sites n = %d (purple)", sum(plot_df$is_spliceosome)))
ggsave("figures_kncc/m5_raw_vs_corrected.pdf", p, width = 6, height = 5.5)
cat("Done.\n")
