# =============================================================================
# 99. WEBS-free clinical analyses (KNCC longitudinal + CPTAC)
#     spliceosome protein score / delta vs PFS & OS; trajectory groups
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2); library(survival); library(survminer) })
set.seed(42)

b <- readRDS("results/kncc/dual_omics_bundle.rds")
pm <- b$prot_samples; cli <- b$prot_annot; delta <- b$prot_delta

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
sp <- intersect(splice, rownames(pm))
cli$sp_score <- colMeans(pm[sp, , drop = FALSE], na.rm = TRUE)[cli$sample]

master <- b$master
cli2 <- cli %>% left_join(master %>% dplyr::select(patient, PFS, OS, `death.0.1.`, AGE, SEX),
                          by = "patient", suffix = c("", ".m")) %>%
  mutate(OS = as.numeric(OS), PFS = as.numeric(PFS),
         death = as.integer(`death.0.1.`), AGE = as.numeric(AGE))

sink("results/reframe/clinical_kncc.txt", split = TRUE)
cat("=== KNCC clinical: spliceosome protein features vs outcomes ===\n\n")

# ---- T1 score vs OS/PFS (primary tumors only) ----
t1 <- cli2 %>% filter(timepoint == "T1", !is.na(OS))
cat("T1 samples with OS:", nrow(t1), "\n")
cx_os <- coxph(Surv(OS, death) ~ sp_score + AGE + SEX, data = t1)
s <- summary(cx_os)
cat("T1 spliceosome score ~ OS (adj age/sex): HR =", round(s$coefficients["sp_score", "exp(coef)"], 2),
    " p =", signif(s$coefficients["sp_score", "Pr(>|z|)"], 3), "\n")
t1p <- t1 %>% filter(!is.na(PFS)) %>% mutate(pfs_event = 1L)
cx_pfs <- coxph(Surv(PFS, pfs_event) ~ sp_score + AGE + SEX, data = t1p)
s2 <- summary(cx_pfs)
cat("T1 spliceosome score ~ PFS: HR =", round(s2$coefficients["sp_score", "exp(coef)"], 2),
    " p =", signif(s2$coefficients["sp_score", "Pr(>|z|)"], 3), "\n")

# ---- delta spliceosome vs OS (landmark at recurrence sampling) ----
d_sp <- colMeans(delta[sp, , drop = FALSE], na.rm = TRUE)
dd <- master %>%
  mutate(OS = as.numeric(OS), PFS = as.numeric(PFS),
         death = as.integer(`death.0.1.`), AGE = as.numeric(AGE)) %>%
  filter(!is.na(OS)) %>%
  mutate(d_sp = d_sp[patient])
dd <- dd[!is.na(dd$d_sp), ]
cat("\npatients with delta + OS:", nrow(dd), "\n")
cx_d <- coxph(Surv(OS, death) ~ d_sp + AGE + SEX, data = dd)
s3 <- summary(cx_d)
cat("delta spliceosome ~ OS: HR =", round(s3$coefficients["d_sp", "exp(coef)"], 2),
    " p =", signif(s3$coefficients["d_sp", "Pr(>|z|)"], 3), "\n")
cat("NOTE: T2-sample analyses are subject to immortal-time/landmark considerations;\n")
cat("      delta features use OS from diagnosis and are exploratory.\n")

# ---- proliferation-trajectory groups vs outcomes ----
prolif <- intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                      "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(delta))
d_pro <- colMeans(delta[prolif, , drop = FALSE], na.rm = TRUE)
ter <- cut(d_pro, breaks = quantile(d_pro, 0:3/3, na.rm = TRUE), include.lowest = TRUE,
           labels = c("Prolif-drop", "Prolif-stable", "Prolif-rise"))
dd$traj <- ter[match(dd$patient, names(d_pro))]
dd$traj <- relevel(dd$traj, ref = "Prolif-stable")
kt <- survdiff(Surv(OS, death) ~ traj, data = dd)
cat("\ntrajectory groups vs OS: log-rank p =", signif(1 - pchisq(kt$chisq, 2), 3), "\n")
print(dd %>% group_by(traj) %>% summarise(n = n(), median_OS = median(OS), .groups = "drop"))

kt2 <- survdiff(Surv(PFS, rep(1L, sum(!is.na(dd$PFS)))) ~ traj,
                data = dd %>% filter(!is.na(PFS)))
cat("trajectory groups vs PFS: log-rank p =",
    signif(1 - pchisq(kt2$chisq, 2), 3), "\n")
sink()

# KM figure
fit_km <- survfit(Surv(OS, death) ~ traj, data = dd)
p <- ggsurvplot(fit_km, data = dd, pval = TRUE, risk.table = TRUE,
                palette = c("#1f77b4", "grey50", "#d62728"),
                legend.title = "Proliferation trajectory",
                title = "KNCC: overall survival by proliferation trajectory")
pdf("figures_reframe/kncc_traj_os.pdf", width = 6.5, height = 6)
print(p)
dev.off()

# ---- CPTAC: spliceosome protein score vs survival (WEBS-free) ----
prot <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct", check.names = FALSE)
pmc <- as.matrix(sapply(prot[, -1], as.numeric)); rownames(pmc) <- prot[[1]]
spc <- intersect(splice, rownames(pmc))
cptac_score <- colMeans(pmc[spc, ], na.rm = TRUE)
clin <- tryCatch(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_CLI.tsi"), error = function(e) NULL)
if (is.null(clin)) clin <- tryCatch(read.delim("E:/KIMI/data/HS_CPTAC_GBM_CLI.tsi"), error = function(e) NULL)
if (!is.null(clin)) {
  cn <- colnames(clin)
  os_col <- cn[grepl("OS|overall|survival", cn, ignore.case = TRUE)]
  st_col <- cn[grepl("status|vital|death", cn, ignore.case = TRUE)]
  id_col <- cn[1]
  cat("\nCPTAC clinical columns:", paste(cn[1:8], collapse = ", "), "\n")
  clin$sample <- clin[[id_col]]
  dfj <- data.frame(sample = names(cptac_score), score = cptac_score) %>%
    inner_join(clin, by = "sample")
  if (length(os_col) >= 1 && length(st_col) >= 1) {
    dfj$os <- as.numeric(dfj[[os_col[1]]])
    dfj$ev <- as.integer(as.logical(dfj[[st_col[1]]]) | dfj[[st_col[1]]] %in% c("1", "Dead", "DEAD", "deceased"))
    ok <- !is.na(dfj$os) & !is.na(dfj$ev) & dfj$os > 0
    if (sum(ok) > 30) {
      cx <- coxph(Surv(os, ev) ~ score, data = dfj[ok, ])
      sc <- summary(cx)
      cat("CPTAC spliceosome protein score ~ OS: HR =",
          round(sc$coefficients[1, "exp(coef)"], 2),
          " p =", signif(sc$coefficients[1, "Pr(>|z|)"], 3),
          " (n =", sum(ok), ")\n")
    }
  }
}
cat("Done.\n")
