# Supplementary Materials

**Patient-resolved proteomic evolution of recurrent glioblastoma: a two-axis structure of proliferative-state transition and inflammatory–metabolic remodeling**

Jia Chen, Junhan Yang, Bing Yang, Tianyi Zhou, Peijun Zhao

---

## Supplementary Methods

### S1. Module definitions and frozen axis loadings

The 23 pre-specified modules (21 MSigDB Hallmark sets with ≥8 quantified proteins in the KNCC proteome, plus two literature-curated sets: 14-marker proliferation core and 10-marker neuronal core) are listed in Table S1 with their protein compositions and frozen PC1–PC3 loadings from the KNCC discovery PCA (correlation matrix, scaled). The PC2 axis is dominated by negative loadings on inflammatory/metabolic-stress modules (IL2–STAT5 −0.318, xenobiotic metabolism −0.313, hypoxia −0.294, apoptosis −0.281, TNFα–NFκB −0.276, ROS −0.258) and positive loadings on proliferative machinery (G2M +0.201, MYC targets +0.177, E2F +0.175, proliferation core +0.170); consequently, a more negative PC2 score indicates stronger inflammatory–metabolic change at recurrence.

### S2. HHS cohort processing and projection

The imputed log2 LFQ-iBAQ matrix (6,977 proteins × 134 samples) from Tatari et al. (2022) was used as published (Gaussian downshift imputation by the original authors). Sample types were parsed from sample-name suffixes (N/PM/RM/PNM/RNM); the 40 matched patients each contributed exactly one primary and one recurrent proteome. Per-patient deltas were computed within cohort only. KNCC module definitions and PC2 loadings were frozen before projection; modules required ≥5 shared proteins between cohorts (per-module coverage 56–100%; PC2-loading-weighted coverage 73%). Robustness: restriction to ≥60%-coverage modules (HHS PC2 p = 3 × 10⁻⁴), exclusion of bottom-tercile abundance (imputation-suspect) proteins (p = 4 × 10⁻⁴), and leave-one-patient-out (40 iterations, all p < 5 × 10⁻⁴). Effect sizes are archived with 10,000× patient bootstrap median CIs (Tables S4–S5).

### S3. GeoMx processing and patient-aware statistics

Raw DCC files (NanoString GeoMx Whole Transcriptome Atlas) were parsed with probe configuration Hsa_WTA_v1.0 (18,556/18,814 RTS IDs mapped to 18,419 genes), Q3-normalized and log2-transformed. ROI metadata (patient, timepoint, morphology stain CD64/GFAP) was taken from the GEO series matrix; 79 ROIs from 11 patients were analyzed. Frozen KNCC signatures were scored per ROI with shared genes. Patients, not ROIs, were the statistical unit: ROI scores were averaged within patient × timepoint × stain before testing (CD64–GFAP contrasts paired within patient, n = 17 patient–timepoint combinations; primary–recurrent comparisons used paired tests where ≥3 paired patients existed, unpaired tests otherwise, reported with n).

### S4. GLASS longitudinal RNA analysis

GLASS TPM matrix and the longitudinal pairing table were used. Gene-level log2(TPM+0.5) frozen-module scores were computed per sample; 221 primary–recurrent pairs matched to the expression matrix were tested with paired Wilcoxon per module and BH-FDR correction (Table S6).

---

## Supplementary Tables

| Table | Contents |
|---|---|
| S1 | 23 pre-specified modules: sources, protein compositions, frozen PC1–PC3 loadings |
| S2 | GeoMx ROI table: 79 ROIs with patient, timepoint, stain, and frozen axis/module/proxy scores |
| S3 | Module-level heterogeneity with bootstrap CIs and technical-variance ratios (KNCC) |
| S4 | HHS module effect sizes (median + bootstrap 95% CI, both cohorts) |
| S5 | HHS frozen PC2 projection effect sizes with bootstrap CIs |
| S6 | GLASS per-module paired tests (221 pairs) |
| S7 | Whole-proteome module adjudication (50 Hallmark modules: raw vs adjusted effects, class) |
| S8 | Matched complex controls (paired deltas + matched-null empirical p) |
| S9 | KSEA-style kinase results (raw and parent-corrected, with LOO sensitivity) |
| S10 | GDSC2 drug-sensitivity screen (271 drugs, brain lines) |
| S11 | DepMap spliceosome-gene dependency screen (126 brain lines) |
| S12 | CPTAC cross-sectional concordance (21 module pairs) |
| S13 | Composition-adjusted adjudication (before/after cell-composition proxies) |
| S14 | Cellular-origin proxy correlations (KNCC primary tumors) |
| S15 | KNCC per-protein paired deltas (11,320 proteins) |
| S16 | Phosphosite-level results (raw and parent-corrected, both correction methods) |

## Supplementary Figures

**Figure S1.** Module overlap audit. (a) Jaccard overlap matrix for the 23 pre-specified modules. (b) Delta correlations before and after removal of shared proteins for top correlated pairs. (c) PCA under strict non-overlapping module definitions (PC1 39.4%, PC2 35.7%). (d) Leave-one-patient-out influence on PC variance shares.

**Figure S2.** Phosphoproteome coverage audit. Distribution of paired observations per phosphosite (median 12 pairs/site; 5,749 sites with ≥40 pairs) and per-patient site coverage.

**Figure S3.** Spliceosome per-protein residual analysis and strict matched-permutation null (10,000 abundance- and variance-matched random sets drawn without replacement from non-spliceosome genes).

**Figure S4.** Baseline-prediction cross-validated error (M0/M1/M2) with permutation p-values for the four pre-specified endpoints.

**Figure S5.** Exploratory survival association by proliferation-trajectory tertiles (OS and PFS; landmark-bias caveat).
