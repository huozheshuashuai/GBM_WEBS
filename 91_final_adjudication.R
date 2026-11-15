# =============================================================================
# 91. Formal competing-hypothesis adjudication (H1 vs H2), KNCC + GSE155434
#     KNCC: centered predictors; intercept = adjusted spliceosome effect;
#     VIF, partial R2, LOO, patient bootstrap, proliferation-matched strata,
#     3 independent proliferation scores (14-marker, Hallmark E2F, G2M)
#     GSE155434: same model at mRNA level (independent cohort, exon counts)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(MASS) })
set.seed(42)

obj <- readRDS("results/kncc/proteome_analysis_ready.rds")
delta <- obj$delta; paired_pts <- obj$paired

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
sp <- intersect(splice, rownames(delta))

# three independent proliferation scores
prolif_core <- intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                           "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(delta))
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
get_hm <- function(term) {
  l <- gmt[grepl(paste0("^", term, "\t"), gmt)]
  if (!length(l)) return(character(0))
  strsplit(l, "\t")[[1]][-(1:2)]
}
e2f <- intersect(get_hm("HALLMARK_E2F_TARGETS"), rownames(delta))
g2m <- intersect(get_hm("HALLMARK_G2M_CHECKPOINT"), rownames(delta))
neuron <- intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                      "GAP43","TUBB3","MAP2"), rownames(delta))
immune <- intersect(c("CD14","AIF1","FCGR3A","CD68","CD163","ITGAM","PTPRC",
                      "CSF1R","TMEM119","CX3CR1"), rownames(delta))

sc <- function(g) colMeans(delta[g, , drop = FALSE], na.rm = TRUE)
d_sp <- sc(sp)

sink("results/kncc/final_adjudication.txt", split = TRUE)
cat("=== FORMAL ADJUDICATION: H1 (spliceosome-specific) vs H2 (cell-state) ===\n\n")

cat("--- Check 1: score composition ---\n")
cat("prolif_core n =", length(prolif_core), "| overlap with spliceosome:", length(intersect(prolif_core, sp)), "\n")
cat("E2F n =", length(e2f), "| overlap:", length(intersect(e2f, sp)), "\n")
cat("G2M n =", length(g2m), "| overlap:", length(intersect(g2m, sp)), "\n")
cat("neuronal n =", length(neuron), "| overlap:", length(intersect(neuron, sp)), "\n")
cat("immune n =", length(immune), "| overlap:", length(intersect(immune, sp)), "\n")

# ---- formal model with UNCENTERED predictors ----
# intercept = expected spliceosome delta for a patient with zero change in
# proliferation / neuronal / immune scores -> the adjusted spliceosome effect
run_model <- function(yc, pc, nc, ic, label) {
  df <- data.frame(y = yc, p = pc, n = nc, i = ic)
  fit <- lm(y ~ p + n + i, data = df)
  cf <- summary(fit)$coefficients
  vif_val <- sapply(c("p","n","i"), function(v) {
    r2 <- summary(lm(df[[v]] ~ . , data = df[, setdiff(c("p","n","i"), v), drop = FALSE]))$r.squared
    1 / (1 - r2)
  })
  # partial R2
  sse_f <- sum(resid(fit)^2)
  pr2 <- sapply(c("p","n","i"), function(v) {
    fr <- lm(as.formula(paste("y ~", paste(setdiff(c("p","n","i"), v), collapse = "+"))), data = df)
    (sum(resid(fr)^2) - sse_f) / sum(resid(fr)^2)
  })
  cat("\n[", label, "]\n")
  cat("  intercept (adjusted spliceosome effect) =", round(cf[1, 1], 4),
      " 95% CI [", round(confint(fit)[1, 1], 4), ",", round(confint(fit)[1, 2], 4), "]",
      " p =", signif(cf[1, 4], 3), "\n")
  cat("  beta_prolif =", round(cf[2, 1], 3), " p =", signif(cf[2, 4], 3),
      " | beta_neuron =", round(cf[3, 1], 3), " p =", signif(cf[3, 4], 3),
      " | beta_immune =", round(cf[4, 1], 3), " p =", signif(cf[4, 4], 3), "\n")
  cat("  VIF:", paste(round(vif_val, 2), collapse = ", "), "\n")
  cat("  partial R2: prolif", round(pr2[1], 3), " neuron", round(pr2[2], 3),
      " immune", round(pr2[3], 3), "\n")
  invisible(list(fit = fit, df = df))
}

d_pro <- sc(prolif_core); d_neu <- sc(neuron); d_imm <- sc(immune)
cat("\n--- Check 2: formal paired-change model (KNCC, n =", length(d_sp), "patients) ---\n")
m1 <- run_model(d_sp, d_pro, d_neu, d_imm, "14-marker proliferation")

m2 <- run_model(d_sp, sc(e2f), d_neu, d_imm, "Hallmark E2F targets")
m3 <- run_model(d_sp, sc(g2m), d_neu, d_imm, "Hallmark G2M checkpoint")

# ---- patient bootstrap of intercept ----
cat("\n--- patient bootstrap (2000x), 14-marker model ---\n")
B <- 2000; n <- length(d_sp)
b0 <- numeric(B)
for (b in seq_len(B)) {
  idx <- sample(n, n, TRUE)
  fb <- lm(d_sp[idx] ~ d_pro[idx] + d_neu[idx] +
             d_imm[idx])
  b0[b] <- coef(fb)[1]
}
ci <- quantile(b0, c(0.025, 0.975), na.rm = TRUE)
cat("intercept bootstrap 95% CI: [", round(ci[1], 4), ",", round(ci[2], 4), "]\n")
cat("raw spliceosome delta (no adjustment): median", round(median(d_sp), 3),
    " | adjusted intercept:", round(coef(m1$fit)[1], 4), "\n")

# ---- leave-one-patient-out ----
loo <- sapply(seq_len(n), function(j) {
  fb <- lm(d_sp[-j] ~ d_pro[-j] + d_neu[-j] +
             d_imm[-j])
  coef(fb)[1]
})
cat("LOO intercept range: [", round(min(loo), 4), ",", round(max(loo), 4), "]\n")

# ---- proliferation-matched strata ----
cat("\n--- proliferation-matched strata (tertiles of delta-proliferation) ---\n")
ter <- cut(d_pro, breaks = quantile(d_pro, 0:3/3), include.lowest = TRUE, labels = c("drop","stable","rise"))
for (lev in levels(ter)) {
  w <- ter == lev
  wt <- suppressWarnings(wilcox.test(d_sp[w], mu = 0))
  cat("  strata", lev, "(n =", sum(w), "): median d_sp =", round(median(d_sp[w]), 3),
      " p =", signif(wt$p.value, 3), "\n")
}
sink()
cat("KNCC part done.\n")

# ================= GSE155434 mRNA-level replication =================
cat("\nbuilding GSE155434 gene counts...\n")
files <- list.files("E:/KIMI/data/GSE155434", pattern = "feature_count_exon.txt.gz$", full.names = TRUE)
sm <- readLines(gzfile("E:/KIMI/data/GSE155434_series_matrix.txt.gz"))
gsm_ids <- strsplit(gsub('"', "", sm[grepl("^!Sample_geo_accession", sm)]), "\t")[[1]][-1]
chars <- sm[grepl("^!Sample_characteristics_ch1", sm)]
pick <- function(pat) sub(paste0(".*", pat, "\\s*"), "", gsub('"', "", strsplit(chars[grepl(pat, chars)], "\t")[[1]][-1]))
smeta <- data.frame(gsm = gsm_ids,
                    ksad = pick("sample id:"),
                    status = pick("primary/recurrent:"),
                    patient = pick("patient id:"),
                    idh = pick("idh mutation:"))
# match count files by GSM prefix (two file flavours: *_feature_count_exon / *_feature_count)
allf <- list.files("E:/KIMI/data/GSE155434", pattern = "feature_count.*\\.txt\\.gz$", full.names = TRUE)
smeta$file <- sapply(smeta$gsm, function(g) {
  hit <- allf[grepl(paste0("^", g, "_"), basename(allf))]
  if (length(hit)) hit[1] else NA_character_
})
smeta <- smeta[!is.na(smeta$file) & smeta$patient != "NA", ]
cat("samples with count files:", nrow(smeta), "\n")

# patients with both Primary and Recurrent (first occurrence each)
pairs <- lapply(split(smeta, smeta$patient), function(d) {
  p <- d[d$status == "Primary", ]; r <- d[d$status == "Recurrent", ]
  if (nrow(p) >= 1 && nrow(r) >= 1) data.frame(patient = d$patient[1], f1 = p$file[1], f2 = r$file[1])
})
pairs <- do.call(rbind, pairs)
cat("paired patients:", nrow(pairs), "\n")

d0 <- read.delim(gzfile(pairs$f1[1]), comment.char = "#", check.names = FALSE)
genes_ord <- unique(sub("\\..*", "", d0$Geneid))
read_gene <- function(f) {
  d <- read.delim(gzfile(f), comment.char = "#", check.names = FALSE)
  tapply(d[[7]], sub("\\..*", "", d$Geneid), sum)[genes_ord]
}
gm1 <- sapply(pairs$f1, read_gene); gm2 <- sapply(pairs$f2, read_gene)
rownames(gm1) <- rownames(gm2) <- genes_ord
suppressMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
gmap <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = genes_ord,
                keytype = "ENSEMBL", columns = "SYMBOL"))
sym <- gmap$SYMBOL[match(genes_ord, gmap$ENSEMBL)]
keep <- !is.na(sym)
gm <- cbind(gm1[keep, ], gm2[keep, ])
sym <- sym[keep]
gm <- rowsum(gm, sym)
lib <- colSums(gm)
lcpm <- t(log2(t(gm + 0.5) / (lib + 1) * 1e6))
sp_g <- intersect(splice, rownames(lcpm))
prol_g <- intersect(prolif_core, rownames(lcpm))
neu_g <- intersect(neuron, rownames(lcpm))
imm_g <- intersect(immune, rownames(lcpm))
s_all <- function(g) colMeans(lcpm[g, , drop = FALSE])
np <- nrow(pairs)
d_sp_m  <- s_all(sp_g)[(np + 1):(2 * np)] - s_all(sp_g)[1:np]
d_pro_m <- s_all(prol_g)[(np + 1):(2 * np)] - s_all(prol_g)[1:np]
d_neu_m <- s_all(neu_g)[(np + 1):(2 * np)] - s_all(neu_g)[1:np]
d_imm_m <- s_all(imm_g)[(np + 1):(2 * np)] - s_all(imm_g)[1:np]

sink("results/kncc/final_adjudication.txt", append = TRUE, split = TRUE)
cat("\n--- Check 3: GSE155434 mRNA-level (independent cohort, n =", np, "pairs) ---\n")
wt <- suppressWarnings(wilcox.test(d_sp_m, mu = 0))
cat("spliceosome mRNA delta: median", round(median(d_sp_m), 3), " p =", signif(wt$p.value, 3), "\n")
wtp <- suppressWarnings(wilcox.test(d_pro_m, mu = 0))
cat("proliferation mRNA delta: median", round(median(d_pro_m), 3), " p =", signif(wtp$p.value, 3), "\n")
run_model(d_sp_m, d_pro_m, d_neu_m, d_imm_m, "GSE155434 mRNA 14-marker model")
sink()
cat("Done.\n")
