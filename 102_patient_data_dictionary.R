# =============================================================================
# 102. 02_PATIENT_DATA_DICTIONARY.csv: full patient-level data dictionary
#      per patient: IDs, institution, clinical, per-layer per-timepoint
#      availability + sample IDs + TMT plex + RNA runs + per-patient missingness
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr) })

b <- readRDS("results/kncc/dual_omics_bundle.rds")
pm <- b$prot_samples; cli <- b$prot_annot
phm <- b$phos_samples; ph_annot <- b$phos_annot
cli_tab <- read.csv("E:/KIMI/data/PDC000514/KNCC_clinical_tableS1.csv")
cli_tab <- cli_tab[!duplicated(cli_tab$PDC_ID) & grepl("^GBM", cli_tab$PDC_ID), ]

# RNA run accessions from ENA metadata
ena <- read.delim("E:/KIMI/ena_prjna.tsv")
rna <- ena[ena$library_strategy == "RNA-Seq", ]
rna$sample_clean <- sub("_RNA$", "", rna$sample_alias)
id2pt <- c(setNames(cli_tab$PDC_ID, cli_tab$SAMPLE_Primary_ID),
           setNames(cli_tab$PDC_ID, cli_tab$SAMPLE_Recur_ID))
rna$PDC_ID <- id2pt[rna$sample_clean]
rna$tp <- sub(".*(T[12])$", "\\1", rna$sample_clean)
rna <- rna[!is.na(rna$PDC_ID), ]
rna_t1 <- tapply(rna$run_accession[rna$tp == "T1"], rna$PDC_ID[rna$tp == "T1"], paste, collapse = ";")
rna_t2 <- tapply(rna$run_accession[rna$tp == "T2"], rna$PDC_ID[rna$tp == "T2"], paste, collapse = ";")

# proteome per-sample annotation (sample id + plex + missingness)
pa <- cli %>% mutate(PDC_ID = sub("KNCC_", "", patient),
                     miss = colMeans(is.na(pm))[sample]) %>%
  dplyr::select(PDC_ID, timepoint, prot_sample = sample, prot_plex = plex, prot_miss = miss)
prot_t1 <- pa %>% filter(timepoint == "T1")
prot_t2 <- pa %>% filter(timepoint == "T2")

pha <- ph_annot %>% mutate(PDC_ID = sub("KNCC_", "", patient),
                           miss = colMeans(is.na(phm))[sample]) %>%
  dplyr::select(PDC_ID, timepoint, phos_sample = sample, phos_miss = miss)
phos_t1 <- pha %>% filter(timepoint == "T1")
phos_t2 <- pha %>% filter(timepoint == "T2")

dd <- cli_tab %>%
  dplyr::select(PDC_ID, Institution, SAMPLE_Primary_ID, SAMPLE_Recur_ID, SEX, AGE,
                PFS, OS, death = `death.0.1.`) %>%
  left_join(prot_t1, by = "PDC_ID") %>% rename(prot_T1_sample = prot_sample,
            prot_T1_plex = prot_plex, prot_T1_miss = prot_miss) %>%
  left_join(prot_t2, by = "PDC_ID") %>% rename(prot_T2_sample = prot_sample,
            prot_T2_plex = prot_plex, prot_T2_miss = prot_miss) %>%
  left_join(phos_t1, by = "PDC_ID") %>% rename(phos_T1_sample = phos_sample, phos_T1_miss = phos_miss) %>%
  left_join(phos_t2, by = "PDC_ID") %>% rename(phos_T2_sample = phos_sample, phos_T2_miss = phos_miss)
dd <- dd[!duplicated(dd$PDC_ID), ]
dd$rna_T1_runs <- rna_t1[dd$PDC_ID]
dd$rna_T2_runs <- rna_t2[dd$PDC_ID]
dd <- dd %>% mutate(
  prot_paired = !is.na(prot_T1_sample) & !is.na(prot_T2_sample),
  phos_paired = !is.na(phos_T1_sample) & !is.na(phos_T2_sample),
  rna_paired = !is.na(rna_T1_runs) & !is.na(rna_T2_runs),
  triple_paired = prot_paired & phos_paired & rna_paired,
  plex_confounded = prot_paired & (prot_T1_plex == prot_T2_plex)
)
# clean empty columns from xlsx parse artifacts
dd <- dd[, !grepl("^timepoint|^NA", colnames(dd))]
dd <- dd[, colSums(!is.na(dd)) > 0]

write.csv(dd, "02_PATIENT_DATA_DICTIONARY.csv", row.names = FALSE)

sink("results/kncc/dictionary_report.txt", split = TRUE)
cat("=== 02_PATIENT_DATA_DICTIONARY ===\n")
cat("patients:", nrow(dd), "\n")
cat("prot paired:", sum(dd$prot_paired), "| phos paired:", sum(dd$phos_paired),
    "| RNA paired:", sum(dd$rna_paired), "| triple:", sum(dd$triple_paired), "\n")
cat("prot paired with T1,T2 in SAME plex (batch confound):", sum(dd$plex_confounded, na.rm = TRUE), "\n")
cat("median prot missingness T1:", round(median(dd$prot_T1_miss, na.rm = TRUE), 3),
    " T2:", round(median(dd$prot_T2_miss, na.rm = TRUE), 3), "\n")
cat("median phos missingness T1:", round(median(dd$phos_T1_miss, na.rm = TRUE), 3),
    " T2:", round(median(dd$phos_T2_miss, na.rm = TRUE), 3), "\n")
sink()
cat("Done.\n")
