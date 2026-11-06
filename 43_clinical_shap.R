# =============================================================================
# 43. MGMT joint stratification + clinical associations + pathway-level SHAP
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(survival); library(survminer); library(dplyr); library(ggplot2) })
set.seed(42)

PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)

# ---- load CGGA with full clinical ----
load_full <- function(tag) {
  expr <- readRDS(paste0("data/raw/CGGA/", tag, "_expr.rds"))
  clin <- readRDS(paste0("data/raw/CGGA/", tag, "_clinical.rds"))
  raw <- read.delim(if (tag == "cgga693") "data/raw/CGGA/CGGA.mRNAseq_693_clinical.20200506.txt"
                    else "data/raw/CGGA/CGGA.mRNAseq_325_clinical.20200506.txt",
                    check.names = FALSE)
  colnames(raw)[1] <- "sample"
  clin <- merge(clin, raw[, c("sample", grep("MGMT", colnames(raw), value = TRUE))],
                by = "sample", all.x = TRUE)
  g <- intersect(names(PAPER_COEFS), rownames(expr))
  rs <- as.vector(t(expr[g, , drop = FALSE]) %*% PAPER_COEFS[g])
  clin$risk <- rs[match(clin$sample, colnames(expr))]
  clin <- clin[!is.na(clin$risk), ]
  clin$risk_group <- factor(ifelse(clin$risk > median(clin$risk), "High", "Low"), levels = c("Low", "High"))
  clin$MGMT <- ifelse(tolower(clin$MGMTp_methylation_status) == "methylated", "Methylated", "Unmethylated")
  clin
}

sink("results/mgmt_clinical_summary.txt", split = TRUE)
for (tag in c("cgga693", "cgga325")) {
  clin <- load_full(tag)
  d <- clin[!is.na(clin$MGMT), ]
  d$joint <- paste(d$MGMT, d$risk_group, sep = " + ")
  cat("\n=== ", tag, " (n =", nrow(d), ") ===\n")
  fit <- survfit(Surv(OS.time, OS) ~ joint, data = d)
  print(surv_pvalue(fit, data = d))
  print(table(d$joint))

  p <- ggsurvplot(fit, data = d, pval = TRUE,
                  title = paste0(toupper(tag), ": MGMT x signature joint stratification"),
                  palette = c("#1f77b4", "#7fb3d5", "#C0392B", "#E8685E"),
                  xlab = "Time (months)", ylab = "Overall survival",
                  legend.title = "Group")$plot
  ggsave(paste0("figures_final/fig_mgmt_joint_", tag, ".pdf"), p, width = 7, height = 6)

  # clinical associations (continuous risk score)
  cat("\n-- Clinical associations (Spearman/t-test) --\n")
  cat("Age: rho =", round(cor(clin$Age, clin$risk, use = "complete.obs", method = "spearman"), 3), "\n")
  for (v in c("Gender", "Grade", "PRS_type", "IDH_mutation", "MGMTp_methylation_status")) {
    if (v %in% colnames(clin)) {
      grp <- as.factor(clin[[v]])
      if (nlevels(grp) >= 2) {
        kt <- kruskal.test(risk ~ grp, data = clin)
        cat(sprintf("  %s: KW p = %.3g\n", v, kt$p.value))
      }
    }
  }
}
sink()
cat(readLines("results/mgmt_clinical_summary.txt"), sep = "\n")

# ---- pathway-level SHAP ----
# aggregate gene SHAP values into Hallmark pathways for the six-gene model
suppressMessages(library(fgsea))
shap <- read.csv("results/shap_values.csv")
gmt <- "data/raw/msigdb_hallmark.gmt"
if (!file.exists(gmt)) {
  download.file("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt",
                gmt, mode = "wb", quiet = TRUE)
}
hallmark <- gmtPathways(gmt)
six_genes <- c("HIF1A", "SIRT2", "MMP2", "ITGA5", "HDAC1", "SPP1")

long <- shap %>%
  tidyr::pivot_longer(cols = all_of(six_genes), names_to = "gene", values_to = "shap")

# pathways containing each signature gene
pw_map <- do.call(rbind, lapply(names(hallmark), function(pw) {
  g <- intersect(hallmark[[pw]], six_genes)
  if (length(g) > 0) data.frame(pathway = pw, gene = g)
}))

if (nrow(pw_map) > 0) {
  pw_shap <- long %>%
    inner_join(pw_map, by = "gene") %>%
    group_by(pathway) %>%
    summarise(mean_abs_shap = mean(abs(shap)), .groups = "drop") %>%
    arrange(desc(mean_abs_shap))
  write.csv(pw_shap, "results/pathway_shap.csv", row.names = FALSE)

  cat("\n=== Pathway-level SHAP (genes -> pathways) ===\n")
  print(head(pw_shap, 10))

  p2 <- ggplot(head(pw_shap, 12), aes(x = reorder(pathway, mean_abs_shap), y = mean_abs_shap)) +
    geom_bar(stat = "identity", fill = "#6B4C9A") + coord_flip() +
    theme_bw(base_size = 11) +
    labs(x = NULL, y = "Mean |SHAP| (signature genes within pathway)",
         title = "Pathway-level SHAP for the six-gene model") +
    scale_x_discrete(labels = function(x) gsub("HALLMARK_", "", x))
  ggsave("figures_final/fig_pathway_shap.pdf", p2, width = 7, height = 5.5)
}
cat("\nDone.\n")
