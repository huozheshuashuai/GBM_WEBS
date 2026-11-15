# =============================================================================
# 101. Data feasibility audit + heterogeneity-vs-technical-error pilot
#      (a) patient-level data dictionary: prot/phos/RNA x T1/T2 availability
#      (b) RNA intersection from ENA metadata (no FASTQ download needed)
#      (c) technical error from duplicate aliquots vs inter-patient heterogeneity
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

b <- readRDS("results/kncc/dual_omics_bundle.rds")
cli_tab <- read.csv("E:/KIMI/data/PDC000514/KNCC_clinical_tableS1.csv")
cli_tab <- cli_tab[!duplicated(cli_tab$PDC_ID) & grepl("^GBM", cli_tab$PDC_ID), ]

# ---------- (a) + (b) patient x layer x timepoint availability ----------
# RNA availability from ENA run metadata
ena <- read.delim("E:/KIMI/ena_prjna.tsv")
rna_meta <- ena[ena$library_strategy == "RNA-Seq", ]
rna_meta$sample_clean <- sub("_RNA$", "", rna_meta$sample_alias)

# map clinical sample IDs <-> patients
id2pt <- setNames(cli_tab$PDC_ID, cli_tab$SAMPLE_Primary_ID)
id2pt <- c(id2pt, setNames(cli_tab$PDC_ID, cli_tab$SAMPLE_Recur_ID))
rna_meta$PDC_ID <- id2pt[rna_meta$sample_clean]
rna_meta$tp <- sub(".*(T[12])$", "\\1", rna_meta$sample_clean)
rna_ok <- rna_meta[!is.na(rna_meta$PDC_ID), ]
cat("RNA runs mapped to patients:", nrow(rna_ok), "/", nrow(rna_meta), "\n")

prot <- b$master
pts <- prot$patient
avail <- data.frame(patient = pts, PDC_ID = prot$PDC_ID,
                    prot_T1 = prot$prot_T1, prot_T2 = prot$prot_T2,
                    phos_T1 = prot$phos_T1, phos_T2 = prot$phos_T2)
rna_avail <- rna_ok %>% group_by(PDC_ID) %>%
  summarise(rna_T1 = any(tp == "T1"), rna_T2 = any(tp == "T2"), .groups = "drop")
avail <- avail %>% left_join(rna_avail, by = "PDC_ID")
avail$rna_T1[is.na(avail$rna_T1)] <- FALSE
avail$rna_T2[is.na(avail$rna_T2)] <- FALSE
avail <- avail %>% left_join(cli_tab %>% dplyr::select(PDC_ID, Institution, SEX, AGE, PFS, OS),
                             by = "PDC_ID")

sink("results/kncc/feasibility_report.txt", split = TRUE)
cat("=== KNCC data feasibility audit ===\n\n")
cat("total patients in master:", nrow(avail), "\n")
cat("proteome paired (T1+T2):        ", sum(avail$prot_T1 & avail$prot_T2), "\n")
cat("phospho paired:                 ", sum(avail$phos_T1 & avail$phos_T2), "\n")
cat("RNA paired (T1+T2):             ", sum(avail$rna_T1 & avail$rna_T2), "\n")
cat("prot+RNA paired:                ", sum(avail$prot_T1 & avail$prot_T2 & avail$rna_T1 & avail$rna_T2), "\n")
cat("prot+phos+RNA paired:           ", sum(avail$prot_T1 & avail$prot_T2 & avail$phos_T1 & avail$phos_T2 &
                                           avail$rna_T1 & avail$rna_T2), "\n")
cat("clinical (PFS+OS) coverage:     ", sum(!is.na(avail$PFS) & !is.na(avail$OS)), "\n")
cat("\nby institution:\n")
print(avail %>% group_by(Institution) %>%
        summarise(n = n(), prot_paired = sum(prot_T1 & prot_T2),
                  phos_paired = sum(phos_T1 & phos_T2),
                  rna_paired = sum(rna_T1 & rna_T2), .groups = "drop"))
sink()
write.csv(avail, "results/kncc/data_dictionary.csv", row.names = FALSE)

# ---------- (c) technical error from duplicate aliquots ----------
# re-read raw proteome keeping duplicate columns separate
f <- "E:/KIMI/data/PDC000514/KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv"
hdr <- strsplit(readLines(f, n = 1), "\t")[[1]]
log_idx <- which(grepl("^Log KNCC_GBM", hdr))
cc <- rep("NULL", length(hdr)); cc[1] <- "character"; cc[log_idx] <- "numeric"
pr <- read.delim(f, sep = "\t", header = TRUE, colClasses = cc, check.names = FALSE)
pr <- pr[!pr$Gene %in% c("Mean", "Median", "StdDev") & !is.na(pr$Gene), ]
cn_all <- colnames(pr)[-1]
cn_base <- sub(":.*", "", sub("^Log ", "", cn_all))
dup_names <- unique(cn_base[duplicated(cn_base)])
cat("\nduplicate-aliquot samples for error estimation:", length(dup_names), "\n")

tech_var <- c()
for (s in dup_names) {
  idx <- which(cn_base == s)
  m <- as.matrix(pr[, -1][, idx, drop = FALSE])
  mode(m) <- "numeric"
  dff <- m[, 1] - m[, 2]
  tech_var <- c(tech_var, var(dff, na.rm = TRUE))
}
mean_tech_var <- mean(tech_var)
cat("per-protein technical variance (within-sample, cross-plex):", round(mean_tech_var, 4), "\n")
cat("  => technical SD ~", round(sqrt(mean_tech_var), 3), "log2 units\n")

# biological heterogeneity: variance of per-patient deltas for predefined modules
delta <- b$prot_delta
source("08_cptac_proteomics.R")
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
get_hm <- function(term) {
  l <- gmt[grepl(paste0("^", term, "\t"), gmt)]
  if (!length(l)) return(character(0))
  intersect(strsplit(l, "\t")[[1]][-(1:2)], rownames(delta))
}
mods <- list(
  Spliceosome = intersect(get_spliceosome_genes(), rownames(delta)),
  E2F = get_hm("HALLMARK_E2F_TARGETS"),
  G2M = get_hm("HALLMARK_G2M_CHECKPOINT"),
  ROS = get_hm("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"),
  IL2_STAT5 = get_hm("HALLMARK_IL2_STAT5_SIGNALING"),
  Neuronal = intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2","GAP43","TUBB3","MAP2"),
                       rownames(delta))
)

sink("results/kncc/feasibility_report.txt", append = TRUE, split = TRUE)
cat("\n=== heterogeneity pilot: is patient-to-patient delta variance >> technical noise? ===\n")
cat("technical variance per protein (from duplicate aliquots):", round(mean_tech_var, 4), "\n")
cat("note: module-score technical variance shrinks by 1/n_genes under independence\n\n")
het <- lapply(names(mods), function(mn) {
  g <- mods[[mn]]
  d <- colMeans(delta[g, , drop = FALSE], na.rm = TRUE)
  n_pt <- sum(!is.na(d))
  v_obs <- var(d, na.rm = TRUE)
  v_tech <- mean_tech_var / length(g)   # technical var of the module mean (upper bound)
  data.frame(module = mn, n_genes = length(g), n_patients = n_pt,
             var_delta = round(v_obs, 4), tech_var_bound = round(v_tech, 5),
             ratio = round(v_obs / v_tech, 1),
             sd_delta = round(sd(d, na.rm = TRUE), 3))
})
het <- do.call(rbind, het)
print(het, row.names = FALSE)
cat("\ninterpretation: ratio >> 1 means patient heterogeneity dominates technical noise\n")
sink()

# figure
het$module <- factor(het$module, levels = het$module[order(-het$ratio)])
p <- ggplot(het, aes(module, ratio)) +
  geom_col(fill = "#7B4EA3") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_text(aes(label = paste0("SD=", sd_delta)), vjust = -0.4, size = 3) +
  theme_bw(base_size = 11) +
  labs(y = "Patient delta variance / technical variance bound",
       title = "Is inter-patient heterogeneity real? (variance decomposition pilot)")
ggsave("figures_kncc/feasibility_heterogeneity.pdf", p, width = 7, height = 4.5)
cat("Done.\n")
