# =============================================================================
# 55. Pseudoreplication control: patient-level pseudobulk + mixed-effects models
#     Cells from the same patient are NOT independent replicates.
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(Seurat); library(lme4); library(lmerTest); library(dplyr); library(ggplot2) })

mal <- readRDS("data/processed/GSE138794_malignant.rds")
md <- mal@meta.data
md$patient <- md$sample   # GSE138794: sample = patient-level library
cat("Malignant cells:", nrow(md), " from", length(unique(md$patient)), "patients\n")

sink("results/pseudobulk_validation.txt", split = TRUE)
cat("=== Pseudoreplication control (patient as unit / random effect) ===\n")
cat("Cells:", nrow(md), " Patients:", length(unique(md$patient)), "\n\n")

# ---------- 1. Neftel state: dominant state per cell ----------
md$Neftel_state <- apply(md[, c("Neftel_MES","Neftel_AC","Neftel_NPC","Neftel_OPC")], 1,
                         function(x) c("MES","AC","NPC","OPC")[which.max(x)])
# pseudobulk: mean WEBS per patient x state (states with >= 20 cells)
pb <- md %>% group_by(patient, Neftel_state) %>%
  summarise(mean_WEBS = mean(WEBS), n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 20)
cat("-- Pseudobulk WEBS by Neftel state (patient x state means) --\n")
print(pb %>% group_by(Neftel_state) %>%
        summarise(n_patients = n(), median_WEBS = round(median(mean_WEBS), 3)))
m1 <- lmer(mean_WEBS ~ Neftel_state + (1 | patient), data = pb)
cat("\nMixed model (pseudobulk, patient as random effect):\n")
print(anova(m1))

# cell-level mixed model for comparison
m1c <- lmer(WEBS ~ Neftel_state + (1 | patient), data = md)
cat("\nCell-level mixed model (patient random effect):\n")
print(anova(m1c))

# ---------- 2. Region type ----------
md_r <- md[md$Region_Type %in% c("Core-like","Margin-like","Normal-like"), ]
md_r$Region_Type <- factor(md_r$Region_Type, levels = c("Normal-like","Core-like","Margin-like"))
cat("\n-- Region type --\n")
pb_r <- md_r %>% group_by(patient, Region_Type) %>%
  summarise(mean_WEBS = mean(WEBS), n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 20)
m2 <- lmer(mean_WEBS ~ Region_Type + (1 | patient), data = pb_r)
cat("Pseudobulk mixed model:\n"); print(anova(m2))
m2c <- lmer(WEBS ~ Region_Type + (1 | patient), data = md_r)
cat("Cell-level mixed model:\n"); print(anova(m2c))

# patient-level paired test: Margin vs Core (patients with both)
pb_w <- pb_r %>% filter(Region_Type %in% c("Core-like","Margin-like")) %>%
  tidyr::pivot_wider(names_from = Region_Type, values_from = mean_WEBS,
                     id_cols = patient, values_fn = mean)
pb_w <- pb_w[complete.cases(pb_w), ]
wt <- wilcox.test(pb_w$`Margin-like`, pb_w$`Core-like`, paired = TRUE)
cat(sprintf("\nPaired patient-level Margin vs Core: n = %d patients, median diff = %.3f, Wilcoxon p = %.3g\n",
            nrow(pb_w), median(pb_w$`Margin-like` - pb_w$`Core-like`), wt$p.value))

# ---------- 3. Hypoxia decoupling within malignant cells ----------
cat("\n-- Hypoxia decoupling (per-patient) --\n")
per_pat <- md %>% group_by(patient) %>%
  summarise(rho = ifelse(n() > 100 && sd(WEBS) > 0 && sd(Hallmark_Hypoxia) > 0,
                         cor(WEBS, Hallmark_Hypoxia, method = "spearman"), NA),
            n_cells = n(), .groups = "drop") %>% filter(!is.na(rho))
wt2 <- wilcox.test(per_pat$rho, mu = 0)
cat(sprintf("Per-patient Spearman rho (WEBS vs Hypoxia): median = %.3f, n = %d patients, one-sample Wilcoxon p = %.3g\n",
            median(per_pat$rho), nrow(per_pat), wt2$p.value))
m3 <- lmer(Hallmark_Hypoxia ~ WEBS + (1 | patient), data = md)
cat("Cell-level mixed model (Hypoxia ~ WEBS + (1|patient)):\n")
print(summary(m3)$coefficients)

# ---------- 4. Pseudotime ----------
pt <- readRDS("data/processed/GSE138794_malignant_pt.rds")
mdp <- pt@meta.data
pt_cols <- grep("pseudotime|sling", colnames(mdp), ignore.case = TRUE, value = TRUE)
cat("\n-- Pseudotime (per-patient) --\n")
cat("pseudotime cols:", paste(pt_cols, collapse = ", "), "\n")
ptcol <- pt_cols[1]
mdp$patient <- mdp$sample
pp <- mdp %>% filter(!is.na(.data[[ptcol]])) %>% group_by(patient) %>%
  summarise(rho = ifelse(n() > 100 && sd(WEBS) > 0 && sd(.data[[ptcol]]) > 0,
                         cor(WEBS, .data[[ptcol]], method = "spearman"), NA),
            .groups = "drop") %>% filter(!is.na(rho))
wt3 <- wilcox.test(pp$rho, mu = 0)
cat(sprintf("Per-patient Spearman rho (WEBS vs pseudotime): median = %.3f, n = %d patients, one-sample Wilcoxon p = %.3g\n",
            median(pp$rho), nrow(pp), wt3$p.value))
m4 <- lmer(WEBS ~ get(ptcol) + (1 | patient), data = mdp[!is.na(mdp[[ptcol]]), ])
cat("Cell-level mixed model (WEBS ~ pseudotime + (1|patient)):\n")
print(summary(m4)$coefficients)
sink()

# ---------- figure: per-patient rho distributions ----------
fig_df <- rbind(
  data.frame(analysis = "Hypoxia decoupling", rho = per_pat$rho),
  data.frame(analysis = "Pseudotime", rho = pp$rho))
p <- ggplot(fig_df, aes(analysis, rho)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_violin(fill = "#EDE7F6", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "#7B4EA3", alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  theme_bw(base_size = 12) +
  labs(title = "Per-patient correlations (pseudoreplication control)",
       subtitle = paste0("Hypoxia: median rho = ", round(median(per_pat$rho), 2),
                         ", p = ", signif(wt2$p.value, 2),
                         " | Pseudotime: median rho = ", round(median(pp$rho), 2),
                         ", p = ", signif(wt3$p.value, 2)),
       x = NULL, y = "Per-patient Spearman rho")
ggsave("figures_final/fig_pseudobulk.pdf", p, width = 7, height = 5)
cat("Done.\n")
