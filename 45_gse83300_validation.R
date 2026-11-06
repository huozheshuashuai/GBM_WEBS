# =============================================================================
# 45. External validation on GSE83300 (50 GBM samples, survival available)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(GEOquery); library(survival); library(survminer); library(dplyr) })

g <- getGEO(filename = "data/raw/GSE83300/GSE83300_series_matrix.txt.gz")
expr <- exprs(g)
pd <- pData(g)

# probe -> symbol via GPL
gpl <- getGEO(g@annotation, destdir = "data/raw/GPL")

fd <- Table(gpl)
fd <- fd[!is.na(fd$ID) & fd$ID != "", ]
rownames(fd) <- fd$ID
sym <- fd[[grep("Gene Symbol|gene_symbol|Symbol", colnames(fd), ignore.case = TRUE)[1]]]
sym <- sub(" ///.*$", "", sym)
common_probes <- intersect(rownames(expr), rownames(fd))
expr <- expr[common_probes, ]; sym <- sym[match(common_probes, rownames(fd))]
keep <- !is.na(sym) & sym != "" & sym != "---"
expr <- expr[keep, ]; sym <- sym[keep]
o <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
expr <- expr[o, ]; sym <- sym[o]
expr <- expr[!duplicated(sym), ]; rownames(expr) <- sym[!duplicated(sym)]
expr <- log2(expr + 1)

# survival
chr_cols <- grep("characteristics", colnames(pd), value = TRUE)
get_field <- function(pat) {
  for (cc in chr_cols) {
    v <- as.character(pd[[cc]])
    if (any(grepl(pat, v, ignore.case = TRUE))) return(sub("^.*?:\\s*", "", v))
  }
  NULL
}
tim <- as.numeric(get_field("overall survival"))
ev <- as.integer(get_field("events"))
clin <- data.frame(sample = rownames(pd), OS.time = tim, OS = ev)
clin <- clin[!is.na(clin$OS.time) & !is.na(clin$OS) & clin$OS.time > 0, ]
expr <- expr[, clin$sample]

cat("GSE83300:", ncol(expr), "samples,", sum(clin$OS), "events\n")

PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)
WRITERS <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP")
ERASERS <- c("HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")

clin$risk <- as.vector(t(expr[intersect(names(PAPER_COEFS), rownames(expr)), ]) %*%
                         PAPER_COEFS[intersect(names(PAPER_COEFS), rownames(expr))])
clin$risk_group <- factor(ifelse(clin$risk > median(clin$risk), "High", "Low"), levels = c("Low", "High"))
clin$WEBS <- colMeans(expr[intersect(WRITERS, rownames(expr)), ]) -
             colMeans(expr[intersect(ERASERS, rownames(expr)), ])

sink("results/gse83300_validation.txt", split = TRUE)
cat("=== GSE83300 external validation (n =", nrow(clin), ") ===\n\n")
cat("-- Six-gene signature --\n")
print(summary(coxph(Surv(OS.time, OS) ~ risk, data = clin))$coefficients)
fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = clin)
print(surv_pvalue(fit, data = clin))
cat("\n-- WEBS (bulk) --\n")
print(summary(coxph(Surv(OS.time, OS) ~ WEBS, data = clin))$coefficients)
sink()
cat(readLines("results/gse83300_validation.txt"), sep = "\n")

p <- ggsurvplot(fit, data = clin, pval = TRUE, title = "GSE83300: six-gene signature",
                palette = c("#1f77b4", "#d62728"), legend.labs = c("Low risk", "High risk"),
                xlab = "Time (months)", ylab = "Overall survival")$plot
ggplot2::ggsave("figures_final/fig_gse83300_km.pdf", p, width = 6, height = 5)
