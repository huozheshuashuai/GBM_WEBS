# Generate supplementary tables for the CL package
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(dplyr))
dir.create("submission_CancerLetters/04_Supplementary/tables", recursive = TRUE, showWarnings = FALSE)

# S1: 23 modules + frozen PC2 loadings
b <- readRDS("results/kncc/dual_omics_bundle.rds")
kdelta <- b$prot_delta
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f2 <- strsplit(l, "\t")[[1]]; f2[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
keep_mods <- c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_DNA_REPAIR",
               "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS",
               "HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
               "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_INTERFERON_GAMMA_RESPONSE",
               "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
               "HALLMARK_ANGIOGENESIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_PI3K_AKT_MTOR_SIGNALING",
               "HALLMARK_P53_PATHWAY","HALLMARK_APOPTOSIS","HALLMARK_HYPOXIA",
               "HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_PROTEIN_SECRETION",
               "HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_MYC_TARGETS_V1")
mods_k <- lapply(hm[keep_mods], function(g) intersect(g, rownames(kdelta)))
mods_k$NEURONAL_CURATED <- intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                                       "GAP43","TUBB3","MAP2"), rownames(kdelta))
mods_k$PROLIF_CURATED <- intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                                     "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(kdelta))
mods_k <- mods_k[sapply(mods_k, length) >= 8]
Md <- t(sapply(mods_k, function(g) colMeans(kdelta[g, , drop = FALSE], na.rm = TRUE)))
Md <- Md[, colSums(is.na(Md)) == 0]
pc <- prcomp(t(Md), scale. = TRUE)
ve <- pc$sdev^2 / sum(pc$sdev^2)

s1 <- data.frame(module = names(mods_k),
                 n_proteins_kncc = sapply(mods_k, length),
                 source = ifelse(grepl("CURATED", names(mods_k)), "curated (literature)", "MSigDB Hallmark"),
                 genes = sapply(mods_k, paste, collapse = ";"),
                 PC1_loading = round(pc$rotation[, 1], 4),
                 PC2_loading = round(pc$rotation[, 2], 4),
                 PC3_loading = round(pc$rotation[, 3], 4))
write.csv(s1, "submission_CancerLetters/04_Supplementary/tables/S1_modules_frozen_loadings.csv",
          row.names = FALSE)

# S2: GeoMx ROI table
gm <- readRDS("results/geomx/geomx_bundle.rds")
roi <- gm$roi_meta %>%
  dplyr::select(gsm, title, patient, timepoint, stain, dcc_stub, PC2,
                mod_IL2STAT5, mod_ROS, mod_EMT, px_Microglia, px_Endothelial) %>%
  rename(ROI_dcc_file = dcc_stub, PC2_axis = PC2, IL2STAT5 = mod_IL2STAT5,
         ROS = mod_ROS, EMT = mod_EMT, Microglia = px_Microglia, Endothelial = px_Endothelial)
write.csv(roi, "submission_CancerLetters/04_Supplementary/tables/S2_geomx_roi_table.csv",
          row.names = FALSE)

# S3: copy key result tables
tabs <- c(
  "results/kncc/heterogeneity_all.csv",
  "results/hhs/hhs_module_effect_sizes.csv",
  "results/hhs/hhs_pc2_effect_sizes.csv",
  "results/kncc/glass_second_axis.csv",
  "results/kncc/module_adjudication.csv",
  "results/kncc/complex_controls.csv",
  "results/kncc/ksea_kinases.csv",
  "results/reframe/gdsc2_spliceosome_drugs.csv",
  "results/reframe/depmap_spliceosome_dependency.csv",
  "results/kncc/cptac_matched_validation.csv",
  "results/kncc/composition_sensitivity.csv",
  "results/kncc/cellular_origin_kncc.csv",
  "results/kncc/paired_differential_proteins.csv",
  "results/kncc/phosphosite_paired.csv"
)
out_names <- c("S3_heterogeneity.csv","S4_hhs_module_effects.csv","S5_hhs_pc2_projection.csv",
               "S6_glass_modules.csv","S7_module_adjudication.csv","S8_complex_controls.csv",
               "S9_ksea_kinases.csv","S10_gdsc2_drugs.csv","S11_depmap_dependency.csv",
               "S12_cptac_concordance.csv","S13_composition_sensitivity.csv",
               "S14_cellular_origin.csv","S15_kncc_perprotein_deltas.csv","S16_phosphosites.csv")
for (i in seq_along(tabs)) {
  if (file.exists(tabs[i])) file.copy(tabs[i], file.path("submission_CancerLetters/04_Supplementary/tables", out_names[i]), overwrite = TRUE)
}
cat("supplementary tables written:", length(list.files("submission_CancerLetters/04_Supplementary/tables")), "\n")
