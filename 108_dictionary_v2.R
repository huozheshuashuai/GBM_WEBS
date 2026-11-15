# =============================================================================
# 108. 02_PATIENT_DATA_DICTIONARY_V2.csv
#      V2 changes: plex_confounded -> same_plex_pair; add mapping evidence,
#      age flags (pediatric GBM0065/GBM0116), access dates, data versions,
#      PFS/OS unit documentation, phospho plex status, fixed analysis sets
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr) })

b <- readRDS("results/kncc/dual_omics_bundle.rds")
pm <- b$prot_samples; cli <- b$prot_annot
phm <- b$phos_samples; ph_annot <- b$phos_annot
cli_tab <- read.csv("E:/KIMI/data/PDC000514/KNCC_clinical_tableS1.csv")
cli_tab <- cli_tab[!duplicated(cli_tab$PDC_ID) & grepl("^GBM", cli_tab$PDC_ID), ]

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

pa <- cli %>% mutate(PDC_ID = sub("KNCC_", "", patient),
                     miss = colMeans(is.na(pm))[sample]) %>%
  dplyr::select(PDC_ID, timepoint, prot_sample = sample, prot_plex = plex, prot_miss = miss)
prot_t1 <- pa %>% filter(timepoint == "T1") %>% dplyr::select(-timepoint)
prot_t2 <- pa %>% filter(timepoint == "T2") %>% dplyr::select(-timepoint)

pha <- ph_annot %>% mutate(PDC_ID = sub("KNCC_", "", patient),
                           miss = colMeans(is.na(phm))[sample]) %>%
  dplyr::select(PDC_ID, timepoint, phos_sample = sample, phos_miss = miss)
phos_t1 <- pha %>% filter(timepoint == "T1") %>% dplyr::select(-timepoint)
phos_t2 <- pha %>% filter(timepoint == "T2") %>% dplyr::select(-timepoint)

dd <- cli_tab %>%
  dplyr::select(PDC_ID, Institution, SAMPLE_Primary_ID, SAMPLE_Recur_ID, SEX, AGE,
                PFS, OS, death = `death.0.1.`) %>%
  mutate(AGE_num = suppressWarnings(as.numeric(AGE)),
         pediatric_flag = !is.na(AGE_num) & AGE_num < 18,
         mapping_evidence = "PDC_ID in proteome/phospho matrix colnames; SAMPLE_*_ID to ENA sample_alias (suffix _RNA stripped); clinical from Kim2024 Table S1",
         PFS_unit = "days", OS_unit = "days",
         IDH_status = NA_character_, MGMT_status = NA_character_,
         treatment = "standard (per CPTAC-3 protocol, see Kim2024)",
         hypermutator = NA_character_) %>%
  left_join(prot_t1 %>% rename_with(~ paste0("prot_T1_", .x), -PDC_ID), by = "PDC_ID") %>%
  left_join(prot_t2 %>% rename_with(~ paste0("prot_T2_", .x), -PDC_ID), by = "PDC_ID") %>%
  left_join(phos_t1 %>% rename_with(~ paste0("phos_T1_", .x), -PDC_ID), by = "PDC_ID") %>%
  left_join(phos_t2 %>% rename_with(~ paste0("phos_T2_", .x), -PDC_ID), by = "PDC_ID")
dd <- dd[!duplicated(dd$PDC_ID), ]
dd$rna_T1_runs <- rna_t1[dd$PDC_ID]
dd$rna_T2_runs <- rna_t2[dd$PDC_ID]
dd <- dd %>% mutate(
  prot_paired = !is.na(prot_T1_prot_sample) & !is.na(prot_T2_prot_sample),
  phos_paired = !is.na(phos_T1_phos_sample) & !is.na(phos_T2_phos_sample),
  rna_paired = !is.na(rna_T1_runs) & !is.na(rna_T2_runs),
  triple_paired = prot_paired & phos_paired & rna_paired,
  same_plex_pair = prot_paired & (prot_T1_prot_plex == prot_T2_prot_plex),
  phos_plex = "not available in processed download (requires PDC aliquot manifest query)",
  sample_region = "unknown (not recorded in Table S1 or PDC processed files)",
  data_version = "PDC000514/PDC000515 current release; downloaded 2026-11-14; clinical Table S1 from PMC10939876 (accessed 2026-11-14); RNA metadata from ENA filereport (accessed 2026-11-14)"
)
dd <- dd[, colSums(!is.na(dd)) > 0]

write.csv(dd, "02_PATIENT_DATA_DICTIONARY_V2.csv", row.names = FALSE)

sink("results/kncc/dictionary_v2_report.txt", split = TRUE)
cat("=== 02_PATIENT_DATA_DICTIONARY_V2 ===\n")
cat("patients:", nrow(dd), "\n")
cat("prot paired:", sum(dd$prot_paired), "| phos paired:", sum(dd$phos_paired),
    "| RNA paired:", sum(dd$rna_paired), "| triple:", sum(dd$triple_paired), "\n")
cat("same-plex pairs (batch-controlled deltas):", sum(dd$same_plex_pair, na.rm = TRUE), "\n")
cat("cross-plex pairs (sensitivity flagged):", sum(dd$prot_paired & !dd$same_plex_pair, na.rm = TRUE), "\n")
cat("pediatric flags (age<18):", sum(dd$pediatric_flag, na.rm = TRUE),
    "->", paste(dd$PDC_ID[dd$pediatric_flag %in% TRUE], collapse = ", "), "\n")
cat("median prot missingness T1/T2:", round(median(dd$prot_T1_prot_miss, na.rm = TRUE), 3), "/",
    round(median(dd$prot_T2_prot_miss, na.rm = TRUE), 3), "\n")
sink()

# field documentation
note_for <- function(fn) {
  dplyr::case_when(
    fn == "PDC_ID" ~ "PDC case ID",
    fn == "Institution" ~ "tissue source site",
    fn == "SAMPLE_Primary_ID" ~ "primary sample ID (institution naming)",
    fn == "SAMPLE_Recur_ID" ~ "recurrent sample ID",
    fn == "SEX" ~ "sex",
    fn == "AGE" ~ "age at diagnosis (years)",
    fn == "PFS" ~ "progression-free survival, days",
    fn == "OS" ~ "overall survival, days",
    fn == "death" ~ "vital status 1=dead",
    fn == "AGE_num" ~ "age numeric",
    fn == "pediatric_flag" ~ "age<18 (GBM0065 age 1, GBM0116 age 17; verify histology, excluded in sensitivity)",
    fn == "mapping_evidence" ~ "ID mapping provenance",
    fn == "PFS_unit" ~ "days",
    fn == "OS_unit" ~ "days",
    fn == "IDH_status" ~ "not in Table S1 (in WES data, not downloaded)",
    fn == "MGMT_status" ~ "not in Table S1",
    fn == "treatment" ~ "CPTAC-3 protocol, see Kim2024",
    fn == "hypermutator" ~ "in WES data, not downloaded",
    grepl("prot_T[12]_prot_sample", fn) ~ "proteome sample ID",
    grepl("prot_T[12]_prot_plex", fn) ~ "TMT plex",
    grepl("prot_T[12]_prot_miss", fn) ~ "proteome missingness fraction",
    grepl("phos_T[12]_phos_sample", fn) ~ "phospho sample ID",
    grepl("phos_T[12]_phos_miss", fn) ~ "phospho missingness fraction",
    grepl("rna_T[12]_runs", fn) ~ "ENA RNA run accessions",
    fn == "prot_paired" ~ "proteome T1+T2 available",
    fn == "phos_paired" ~ "phospho T1+T2 available",
    fn == "rna_paired" ~ "RNA T1+T2 available",
    fn == "triple_paired" ~ "all three layers paired",
    fn == "same_plex_pair" ~ "T1 and T2 in same TMT plex (deltas batch-controlled)",
    fn == "phos_plex" ~ "phospho plex unavailable in processed files",
    fn == "sample_region" ~ "tumor region not recorded",
    fn == "data_version" ~ "data versions and access dates",
    TRUE ~ "")
}
fields <- data.frame(field = colnames(dd), note = note_for(colnames(dd)))
write.csv(fields, "02_PATIENT_DATA_DICTIONARY_V2_fields.csv", row.names = FALSE)
cat("field notes written.\n")
cat("Done.\n")
