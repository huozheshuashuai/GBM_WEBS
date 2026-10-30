# =============================================================================
# 26. Key reviewer-facing figures
# =============================================================================
# F1: CopyKAT purity (aneuploid fraction per sample + WEBS in aneuploid cells)
# F2: CMap-style drug bar plot (top WEBS-sensitizing drugs, GDSC2)
# F3: DepMap scatter (WEBS vs Chronos: MED1, SF3B1, MLLT1, KAT2A)
# F4: GDSC2 scatter (WEBS vs LN_IC50: Vorinostat, Entinostat)
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(patchwork); library(ggrepel)
})
set.seed(42)
dir.create("figures", showWarnings = FALSE)

# ---------- F1: CopyKAT ----------
obj <- readRDS("data/processed/GSE138794_final.rds")
pred <- readRDS("data/processed/copykat_calls.rds")
rn <- sub("^[^.]+\\.", "", rownames(pred))
rn <- sub("\\.([0-9]+)$", "-\\1", rn)
pred$cell <- rn
obj$copykat <- pred$copykat.pred[match(colnames(obj), pred$cell)]
obj$copykat[obj$copykat == "not.defined"] <- NA

per_sample <- obj@meta.data %>%
  filter(!is.na(copykat)) %>%
  group_by(sample) %>%
  summarise(n_called = n(),
            frac_aneuploid = mean(copykat == "aneuploid"), .groups = "drop")

p1a <- ggplot(per_sample, aes(x = reorder(sample, -frac_aneuploid), y = frac_aneuploid)) +
  geom_bar(stat = "identity", fill = "#d62728", alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Aneuploid fraction (CopyKAT)",
       title = "A. Tumor-cell fraction per sample")

aneu_meta <- obj@meta.data %>% filter(copykat == "aneuploid")
p1b <- ggplot(aneu_meta, aes(x = WEBS_group, y = WEBS, fill = WEBS_group)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.12, outlier.size = 0.3) +
  scale_fill_manual(values = c("LLS" = "#1f77b4", "LDTS" = "#ff7f0e", "LHS" = "#d62728")) +
  theme_bw(base_size = 11) + theme(legend.position = "none") +
  labs(x = NULL, y = "WEBS",
       title = sprintf("B. WEBS in aneuploid-only cells (KW p < 1e-16, n = %d)", nrow(aneu_meta)))

pdf("figures/fig_copykat_purity.pdf", width = 12, height = 5)
print(p1a + p1b)
dev.off()

# ---------- F2: CMap-style drug bar ----------
drug_cor <- read.csv("results/gdsc2_drug_webs_corr.csv")
top_drugs <- drug_cor %>% filter(q < 0.1) %>% arrange(spearman_r) %>% head(20) %>%
  mutate(direction = ifelse(spearman_r < 0, "Sensitizes high-WEBS", "Resistant in high-WEBS"))

p2 <- ggplot(top_drugs, aes(x = reorder(DRUG_NAME, spearman_r), y = spearman_r,
                            fill = direction)) +
  geom_bar(stat = "identity") + coord_flip() +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("Sensitizes high-WEBS" = "#d62728",
                               "Resistant in high-WEBS" = "#1f77b4")) +
  theme_bw(base_size = 11) +
  labs(x = NULL, y = "Spearman rho (WEBS vs LN_IC50)",
       fill = NULL, title = "Top WEBS-sensitizing compounds (GDSC2)")

pdf("figures/fig_cmap_drugs.pdf", width = 7, height = 6)
print(p2)
dev.off()

# ---------- F3: DepMap scatter ----------
webs_cl <- read.csv("results/depmap_cellline_webs.csv")
models <- read.csv("data/raw/DepMap/Model.csv", stringsAsFactors = FALSE)
gbm_lines <- models$ModelID[grepl("glioblastoma|glioma|astrocytoma|gliosarcoma",
                                  models$OncotreePrimaryDisease, ignore.case = TRUE) |
                            grepl("glioblastoma|gliosarcoma", models$OncotreeSubtype,
                                  ignore.case = TRUE)]

# read only needed gene columns from the big CRISPR file
hdr <- read.csv("data/raw/DepMap/CRISPRGeneEffect.csv", nrows = 1, check.names = FALSE)
keep_genes <- c("MED1", "MLLT1", "SF3B1", "KAT2A")
keep_cols <- which(sub(" .*$", "", colnames(hdr)) %in% keep_genes)
cc <- rep("NULL", ncol(hdr)); cc[1] <- "character"; cc[keep_cols] <- "numeric"
crispr <- read.csv("data/raw/DepMap/CRISPRGeneEffect.csv", colClasses = cc,
                   check.names = FALSE)
colnames(crispr) <- c("ModelID", keep_genes[match(colnames(crispr)[-1], colnames(hdr)[keep_cols])])
colnames(crispr)[-1] <- sub(" .*$", "", colnames(hdr)[keep_cols])

sc_df <- webs_cl %>% filter(sample %in% gbm_lines) %>%
  inner_join(crispr, by = c("sample" = "ModelID")) %>%
  tidyr::pivot_longer(cols = all_of(keep_genes), names_to = "gene", values_to = "chronos")

p3 <- ggplot(sc_df, aes(x = WEBS, y = chronos)) +
  geom_point(size = 2, alpha = 0.6, color = "#555555") +
  geom_smooth(method = "lm", se = TRUE, color = "#d62728") +
  geom_hline(yintercept = -1, linetype = "dashed", color = "grey50") +
  facet_wrap(~ gene, ncol = 2, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(y = "Chronos dependency score (more negative = more dependent)",
       title = "WEBS vs CRISPR dependency in GBM cell lines (DepMap 24Q2)")

pdf("figures/fig_depmap_scatter.pdf", width = 9, height = 7)
print(p3)
dev.off()

# ---------- F4: GDSC2 scatter ----------
gdsc <- read.csv("data/raw/GDSC2/GDSC2_fitted_dose_response_24Jul22.csv",
                 stringsAsFactors = FALSE)
m2 <- models[, c("ModelID", "CellLineName")]
gd <- gdsc %>% filter(DRUG_NAME %in% c("Vorinostat", "Entinostat")) %>%
  inner_join(m2, by = c("CELL_LINE_NAME" = "CellLineName")) %>%
  mutate(WEBS = webs_cl$WEBS[match(ModelID, webs_cl$sample)]) %>%
  filter(!is.na(WEBS))

p4 <- ggplot(gd, aes(x = WEBS, y = LN_IC50)) +
  geom_point(size = 1.6, alpha = 0.4, color = "#555555") +
  geom_smooth(method = "lm", se = TRUE, color = "#d62728") +
  facet_wrap(~ DRUG_NAME) +
  theme_bw(base_size = 12) +
  labs(y = "LN IC50 (lower = more sensitive)",
       title = "HDAC inhibitor sensitivity vs WEBS (GDSC2)",
       caption = "Vorinostat: rho = -0.16, q = 1.1e-4 | Entinostat: rho = -0.18, q = 8.7e-5")

pdf("figures/fig_gdsc2_scatter.pdf", width = 10, height = 5)
print(p4)
dev.off()

message("All key figures saved to figures/")
