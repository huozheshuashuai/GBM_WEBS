# Reproducibility re-run battery: core result-generating scripts
# Run in order; each overwrites its results/ outputs
setwd("E:/KIMI/GBM_WEBS_Project")

scripts <- c(
  "08_cptac_proteomics.R",   # CPTAC core: 107/111, rho 0.79
  "20_acetyl_sites.R",       # acetylation OR=7.36
  "29_paper_sixgene_validation.R",  # signature validation
  "28b_specificity_10k.R",   # 10k specificity control
  "44b_virtual_ko_abs.R",    # virtual KO absolute
  "45_gse83300_validation.R",
  "47_glass_validation.R",
  "48_splicing_webs.R",
  "49_phosphoproteome.R",
  "57_ptm_abundance_correction.R",
  "58_uncoupling_rigor.R",
  "59_incremental_value.R",
  "60_splicing_lgg_replication.R",
  "61_signature_performance.R",
  "62_threshold_sensitivity.R",
  "65_glass_firstpair.R"
)

log_con <- file("results/rerun_verification.log", open = "wt")
for (s in scripts) {
  t0 <- Sys.time()
  msg <- tryCatch({
    source(s, local = new.env())
    "OK"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%s  %s  (%.1f min)\n", s, msg,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))), file = log_con)
  flush(log_con)
  cat(s, "->", msg, "\n")
}
close(log_con)
cat("BATTERY DONE\n")
