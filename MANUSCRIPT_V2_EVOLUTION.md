# Title

Patient-resolved proteomic evolution of recurrent glioblastoma: a two-axis structure of proliferative-state transition and inflammatory–metabolic remodeling

# Running title

Patient-specific proteome evolution in recurrent GBM

# Authors

Jia Chen*, Junhan Yang, Bing Yang, Tianyi Zhou, Peijun Zhao

Department of Bioinformatics, Chengdu Neusoft University, Chengdu, China

\* Corresponding author: Jia Chen. Email: 1176625889@qq.com. ORCID: 0009-0002-0082-3050

# Abstract

**Background.** Glioblastoma (GBM) recurs through a molecular evolution that has been characterized at the population level, but whether patients recur through shared or divergent molecular paths — and whether recurrence changes are intrinsic regulatory events or consequences of cell-state transitions — remains unresolved.

**Methods.** We re-analyzed the longitudinal KNCC proteogenomic cohort (105 paired primary–recurrence proteomes; 89 patients with paired phosphoproteomes; clinical annotation for 109 patients) with a patient-resolved, error-aware framework. Technical variance was estimated from cross-plex duplicate aliquots. Paired module-level deltas were decomposed against proliferation, neuronal and immune state scores using uncentered competing-hypothesis models (the intercept estimates state-independent change), with three independent proliferation definitions, patient bootstrap and leave-one-patient-out robustness. Phosphosite changes were tested against parent-protein changes using regression residuals; kinase-substrate enrichment used OmniPath-curated sets with leave-one-substrate-out sensitivity. Cross-sectional support came from CPTAC-GBM (n = 99, paired RNA–protein–PTM).

**Results.** Inter-patient variance of module-level recurrence changes exceeded technical noise by 20- to 346-fold and survived a full sensitivity battery (same-plex pairs, pediatric exclusion, coverage thresholds, leave-one-institution-out). Recurrence evolution was structured along two largely independent axes: a proliferative-to-neuronal axis (PC1, 41% of variance; tightly co-evolving E2F/G2M/MYC/DNA-repair modules, pairwise ρ > 0.9) and an inflammatory–metabolic axis (PC2, 38%; IL2–STAT5, xenobiotic metabolism, hypoxia, ROS). Patients distributed along a continuum without stable discrete subtypes (silhouette < 0.33); hypoxia and EMT modules showed near-50% direction consistency, meaning their null average change masks bidirectional, patient-divergent evolution. The KNCC module co-variation structure replicated in the independent CPTAC-GBM proteome (21 module pairs, ρ = 0.94). Core machinery complexes (spliceosome, ribosomes, RNA polymerase II, DNA replication) contracted concordantly with the proliferation axis, whereas the proteasome was spared (matched-null empirical p = 0.47 vs ≤ 0.001 for machinery complexes). Formal adjudication attributed most module changes to the state transition; after additional cell-composition adjustment, three programs — ROS, IL2–STAT5, and EMT — retained small independent up-regulation. Phosphosite and kinase-level changes (ATM, CDK2, CDK7 down) were entirely accounted for by parent-protein abundance after correction (0/74 kinases FDR < 0.05), and baseline molecular states did not improve prediction of recurrence endpoints beyond clinical covariates and the target's own baseline in nested cross-validation. Proliferation-trajectory subgroups were associated with overall survival (log-rank p = 0.004; exploratory, subject to landmark bias).

**Conclusions.** GBM recurrence follows a two-axis, patient-resolved evolutionary structure in which a tightly coupled proliferative-machinery contraction coexists with a partially independent inflammatory–metabolic program. Most longitudinal proteome and phosphoproteome changes are state-coupled or abundance-driven; competing-hypothesis adjudication is required before attributing module-level changes to independent regulation.

**Keywords:** glioblastoma; recurrence; longitudinal proteomics; phosphoproteomics; patient heterogeneity; cell-state transition; competing hypotheses

# 1. Introduction

Glioblastoma recurs after therapy in essentially all patients. The reference longitudinal proteogenomic study (Kim et al., *Cancer Cell* 2024) established that recurrent GBM shifts toward neuronal-like programs, enriches for specific multi-omic subtypes, and acquires targetable features such as BRAF signaling. That study, like most longitudinal analyses, characterized the *average* trajectory: what changes, on average, between primary and recurrent tumors.

Three questions remain open. First, do all patients traverse the same trajectory, or does recurrence evolution vary across patients beyond technical noise? Second, which recurrence-associated changes are intrinsic regulatory events, and which are statistical shadows of the proliferative-to-neuronal state transition or of protein-abundance change? Third, does the primary tumor's molecular state carry information about the recurrent state beyond its own baseline?

We re-analyze the KNCC longitudinal cohort (105 paired proteomes; 89 paired phosphoproteomes) with a framework designed to answer these questions: technical-error-aware variance decomposition, uncentered competing-hypothesis models that separate state-coupled from state-independent change, matched-complex controls, regression-based parent-protein correction for phosphosites, and nested cross-validated baseline-prediction models with pre-registered stop rules. We report positive and negative results explicitly.

# 2. Methods

## 2.1 Data

KNCC GBM-evolution study: PDC000514 (TMT11 proteome; after QC: 11,320 proteins × 214 tumor samples; 109 patients; 105 paired T1/T2) and PDC000515 (phosphoproteome; 24,015 sites × 185 samples; 96 patients; 89 paired). Clinical variables from the published Table S1 (123 patients; PFS, OS, vital status, age, sex, institution). Cross-sectional support: CPTAC-GBM (LinkedOmics; proteome, acetylome, phosphoproteome, RNA-seq; n = 99). External RNA: GSE155434 (exon-level counts; 11 primary–recurrent pairs). Complex definitions: ComplexPortal (9606) with traceable name-based assignment and overlapping membership allowed. Kinase–substrate sets: OmniPath (PhosphoSite, phosphoELM, PhosphoNetworks, KEGG sources; 4,461 directed edges, 327 kinases, 2,342 substrates).

## 2.2 Quality control and pairing

Non-tumor ("Other"/reference) channels were excluded; duplicate aliquots across plexes were averaged and separately used to estimate technical variance. Of 105 paired patients, 100 had T1 and T2 within the same TMT plex, so paired deltas are batch-controlled by design. Median per-sample missingness: proteome 7.4%; phosphosite pairs median 12 pairs/site (sites with ≥40 pairs: 5,749). Patient-level data dictionary: `02_PATIENT_DATA_DICTIONARY.csv`.

## 2.3 Variance decomposition and heterogeneity

Per-protein technical variance was estimated from seven cross-plex duplicate aliquots. Observed inter-patient variance of per-patient module deltas was compared with the technical-variance bound for the module mean (technical variance ÷ module size; conservative under member correlation). Ratios ≫ 1 indicate biological heterogeneity dominating technical noise.

## 2.4 Module-level longitudinal analysis and competing-hypothesis adjudication

Pre-defined modules (22 Hallmark sets with ≥8 quantified proteins plus curated proliferation and neuronal sets) were scored per sample. Per-patient deltas were tested with paired Wilcoxon. Competing-hypothesis models regressed module deltas on proliferation-, neuronal- and immune-score deltas with uncentered predictors: the intercept estimates the expected module change in patients without state-score change. Robustness: three independent proliferation definitions (14-marker core; Hallmark E2F; G2M), VIF, partial R², patient bootstrap (2,000×), leave-one-patient-out, and proliferation-matched tertile analysis. Marker-set overlap with target modules was excluded.

## 2.5 Co-evolution, trajectory structure, and composition sensitivity

Module-delta correlations used per-patient Spearman with joint patient bootstrap CIs. Trajectory structure was assessed by PCA of the module-delta matrix and by k-means silhouette for k = 2–4, with the continuum model preferred unless discrete clustering was stable. The competing-hypothesis model was extended with cell-composition proxies (Neftel malignant-state modules MES/AC/OPC/NPC and microglia, endothelial, pericyte, oligodendrocyte and neuron marker sets), excluding any genes overlapping the target module. Heterogeneity robustness was assessed across: same-plex pairs (n = 100), pediatric exclusion, ≥80% protein coverage, and leave-one-institution-out analysis.

## 2.6 Matched complex controls

Spliceosome (KEGG), 60S/40S ribosomes, 26S/20S proteasome, RNA polymerase II, DNA-replication (MCM/RF-C) and COPI/COPII complexes were compared by per-protein paired deltas and by empirical p-values against 1,000 abundance-decile-matched random sets drawn without replacement from non-member genes.

## 2.7 Phosphoproteome

Per-site paired deltas (≥40 pairs required) were tested by Wilcoxon. Parent-protein correction used per-site regression of site delta on parent-protein delta (residuals), with delta-subtraction as sensitivity. Gene-level corrected signals were aggregated to Hallmark modules and to kinase-substrate sets (kinases with ≥5 covered substrates; Mann–Whitney of substrates vs non-substrates; BH-FDR; leave-one-substrate-out for top hits).

## 2.8 Baseline prediction

Pre-specified T2 endpoints (neuronal, ROS, IL2–STAT5, proliferation module scores) were modeled as Y_T2 ~ clinical (age, sex) (M0), + target T1 (M1), + other pre-specified T1 module scores (M2). Repeated 10-fold patient-level cross-validation (10 repeats) estimated MSE; M2-vs-M1 gains were tested by patient-level permutation (500×). No feature selection was performed inside folds (fixed predictor sets).

## 2.9 External and supporting analyses

CPTAC-GBM cross-sectional analyses: per-protein protein~mRNA residuals; spliceosome discordance-score spread vs 10,000 abundance- and variance-matched random sets; site-level PTM correlations with the discordance axis, raw and parent-corrected; and cross-sectional replication testing of the KNCC module co-variation structure (21 module-pair Spearman correlations). GSE155434 gene-level counts (summed from exon counts) supported transcript-level state-score comparisons. Proliferation-trajectory groups (tertiles of proliferation-score delta) were compared for OS/PFS by log-rank tests with landmark-bias caveats.

## 2.10 Statistics and reproducibility

R 4.5.2. Patients (not samples or cells) were the statistical unit. Permutation controls excluded the tested set from the background, matched abundance and variance, and sampled without replacement. Multiple testing: Benjamini–Hochberg. Code: https://github.com/huozheshuashuai/GBM_WEBS (scripts 86–107). Generative AI tools were used for language editing and code-formatting assistance; all analyses were executed and verified by the authors.

# 3. Results

## 3.1 A patient-resolved longitudinal proteome resource with quantified technical noise

After QC, 105 patients contributed paired primary–recurrent proteomes (11,320 proteins; median per-sample missingness 7.4%); 89 of them also contributed paired phosphoproteomes, and 70 had additionally paired RNA-seq available in SRA. Three outlier samples (median sample correlation below the 1st percentile) and two pediatric cases (ages 1 and 17) were flagged and handled in sensitivity analyses. Technical variance from cross-plex duplicate aliquots was 0.066 log2² per protein (SD ≈ 0.26). Patient-to-patient variance of module-level recurrence deltas exceeded the technical bound by 20- to 346-fold across the seven frozen modules, and this heterogeneity was unchanged across the full sensitivity battery: same-plex pairs only (n = 100), exclusion of pediatric cases, proteins with ≥80% quantification, and leave-one-institution-out analysis (five institutions) (Fig. 1). Two modules — hypoxia and EMT — showed near-50% direction consistency, i.e., their near-zero mean recurrence change masks bidirectional, patient-divergent evolution; proliferation and neuronal modules were directionally consistent (71% down, 74% up). Inter-patient differences in recurrence evolution are therefore biological, not technical.

## 3.2 Recurrence evolution follows two largely independent axes

Per-patient deltas of 23 pre-defined modules formed a strongly structured co-evolution map (Fig. 2). The proliferative machinery modules (E2F targets, G2M checkpoint, MYC targets, DNA repair) co-varied tightly (pairwise ρ = 0.92–0.96; bootstrap 95% CIs within ±0.05), as did EMT with angiogenesis (ρ = 0.91), glycolysis with hypoxia (ρ = 0.89), IL2–STAT5 with xenobiotic metabolism (ρ = 0.87), and oxidative phosphorylation with neuronal programs (ρ = 0.84). PCA of the module-delta space resolved two dominant axes: PC1 (41% of variance), aligned with the proliferation delta (ρ = 0.72) and anti-aligned with neuronal programs; and PC2 (38%), dominated by inflammatory and metabolic-stress modules (IL2–STAT5, xenobiotic metabolism, hypoxia, apoptosis, ROS) and largely orthogonal to the proliferation axis. Patients occupied a continuum without evidence for stable discrete recurrence subtypes (mean silhouette 0.32/0.28/0.25 for k = 2/3/4).

## 3.3 Most module changes are coupled to the state transition; a small metabolic–inflammatory set is not

Uncentered competing-hypothesis models (module delta ~ proliferation + neuronal + immune deltas) attributed the majority of recurrence change to the state transition: 20 of 50 Hallmark modules with significant raw change showed no independent intercept after adjustment (class A; including E2F, G2M, MYC, DNA repair, unfolded-protein response, OXPHOS). Eleven modules retained small but significant independent intercepts (class B), all positive, dominated by metabolic and inflammatory programs. Adding cell-composition proxies (Neftel malignant states and microglia/endothelial/pericyte/oligodendrocyte/neuron marker sets) to the model removed the hypoxia signal (p = 0.23 after composition adjustment) but preserved the reactive-oxygen-species (p = 5.6 × 10⁻⁴), IL2–STAT5 (p = 5.4 × 10⁻⁴), and EMT (p = 0.003) residuals, indicating that these three programs' recurrence up-regulation is robust to both state-transition and cell-composition explanations (Fig. 3).

## 3.4 The spliceosome as a worked example of adjudication

The spliceosome — a priori candidate for specific dysregulation — declined at recurrence (105/110 proteins FDR < 0.05, direction-uniform; strict matched permutation empirical p = 10⁻⁴). However, the decline was uniform across ComplexPortal subcomplexes (U1, U2, U4/U6, U5, tri-snRNP, NineTeen/Prp19, SR proteins, hnRNPs, DDX/DHX helicases; Kruskal–Wallis p = 0.29), and 80–85% of it was explained by the state transition (adjusted intercept −0.03 to −0.04 log2 across three proliferation definitions; bootstrap CI excluding zero; LOO-stable). Matched complex controls showed the decline belongs to a generic contraction of proliferation-linked machinery (DNA-replication complex −0.22, spliceosome −0.20, 40S ribosome −0.15, RNA polymerase II −0.14, 60S ribosome −0.12; all matched-null p ≤ 0.001), while the proteasome was spared (−0.04; p = 0.47). In independent paired RNA-seq (GSE155434, 11 pairs), the spliceosome showed no transcript-level decline (p = 0.77), consistent with protein-level regulation, though underpowered. We therefore classify the spliceosome decline as predominantly state-coupled with a small unexplained residual of uncertain significance (Fig. 4).

## 3.5 Phosphoproteome changes are abundance shadows after correction

Raw phosphosite analysis found strong recurrence changes, including spliceosome sites (47/107 FDR < 0.05, all down) and kinase-substrate signatures (ATM, CDK2, CDK7 down; CSNK2A1, PRKACA up; raw KSEA FDR < 0.05). After per-site regression on parent-protein deltas, no spliceosome site (0/107), no Hallmark module, and no kinase (0/74) retained FDR-significant independent change; nominally significant kinases were single-substrate-driven on leave-one-substrate-out analysis (Fig. 5). The longitudinal phosphoproteome at this cohort's effective coverage (median 12 pairs/site) therefore provides no evidence for abundance-independent signaling rewiring; raw kinase changes are stoichiometric consequences of protein-abundance change.

## 3.6 Baseline state does not predict recurrence endpoints beyond baseline

Across four pre-specified endpoint modules (neuronal, ROS, IL2–STAT5, proliferation), M2 models (clinical + target baseline + other T1 module scores) did not reduce cross-validated error relative to M1 (clinical + target baseline) — gains were negative (−7% to −14%) with permutation p = 0.42–0.87 (Fig. 6a). The primary tumor's broader molecular state carried no stable incremental information about these recurrence endpoints in this cohort; we therefore make no predictive claim.

## 3.7 External structural replication, exploratory clinical association, and boundary results

The module co-variation structure discovered in KNCC longitudinal deltas replicated strongly in the independent CPTAC-GBM proteome (99 primary tumors): the 21 module-pair correlations matched with ρ = 0.94 (p = 5 × 10⁻⁶), including the proliferation–neuronal anticorrelation (KNCC ρ = −0.39, CPTAC ρ = −0.31), IL2–STAT5–ROS co-variation (0.67 vs 0.78), and hypoxia–EMT co-variation (0.77 vs 0.78) (Fig. 6a). Proliferation-trajectory tertiles (drop/stable/rise) differed in overall survival (median 870/730/945 days; log-rank p = 0.004) and progression-free survival (p = 0.044); because trajectory assignment requires survival to recurrence, these are exploratory and subject to immortal-time/landmark bias (Fig. 6b). CPTAC cross-sectional analysis further confirmed that spliceosome protein–mRNA discordance spread exceeds matched random complexes (empirical p < 10⁻⁴) and that site-level PTM enrichment on the discordance axis inverts after parent-protein correction (acetylation OR 9.7 → 0.43; phosphorylation OR 9.7 → 0.19), independently corroborating the stoichiometric interpretation. Dependency (DepMap) and drug-sensitivity (GDSC2) screens on the spliceosome axis were negative at FDR thresholds and are reported as boundary results.

# 4. Discussion

## 4.1 Principal findings

Re-analyzing 105 paired GBM proteomes with error-aware, patient-resolved statistics, we find that recurrence evolution is (i) genuinely heterogeneous across patients — inter-patient module-delta variance exceeds technical noise by one to two orders of magnitude; (ii) structured along two largely independent axes — a tightly coupled proliferative-to-neuronal axis and a partially independent inflammatory–metabolic axis; (iii) predominantly state-coupled at the module level, with a small set of metabolic/inflammatory modules retaining independent up-regulation; (iv) dominated by protein-abundance change rather than independent signaling rewiring at the phosphosite level; and (v) not predictable, at module resolution, from the primary tumor's broader molecular state beyond its own baseline.

## 4.2 Relationship to the reference study

Kim et al. established the population-average trajectory (neuronal shift, subtype transitions, BRAF) and performed population-level multi-omics comparisons, including mRNA–protein concordance and kinase-network inference (SPHINKS). Our contribution is orthogonal: we move from population averages to patient-resolved evolution with explicit technical-error accounting, formal competing-hypothesis adjudication of state coupling, matched-complex controls, corrected phosphosite analysis with coverage audits, and pre-registered predictive tests with negative results reported. Several population-level findings (e.g., raw kinase changes) are reproduced but shown to be abundance-driven.

## 4.3 The two-axis model

The proliferative axis recapitulates known biology; the second axis — inflammatory and metabolic-stress programs (IL2–STAT5, xenobiotic metabolism, hypoxia, ROS) varying across patients independently of proliferation change — is, to our knowledge, not described in the reference study. Its independence from the state transition and its uniform up-regulation among class-B modules suggest a distinct adaptive layer of recurrence biology. Whether it reflects tumor-intrinsic stress adaptation or microenvironmental remodeling cannot be resolved from bulk proteomes alone; both interpretations are biologically meaningful and are stated as alternatives.

## 4.4 Methodological implications

Three practices made the difference between a naive and an adjudicated reading of this cohort: uncentered state-covariate models (residual-median tests on intercept models are uninformative); matched permutation controls that exclude the tested set, sample without replacement, and match variance as well as abundance; and parent-protein regression correction before any PTM claim. We recommend these as minimal requirements for module-specific longitudinal claims in TMT proteogenomics.

## 4.5 Limitations

This is a secondary analysis of a single published cohort; internal cross-validation is not independent validation, and no external longitudinal proteome cohort of comparable size exists to our knowledge. Phosphoproteome sparsity (median 12 pairs/site) limits site-level power; null corrected results are interpreted as absence of detectable signal, not proof of absence. Module scores compress biology and can mask protein-specific effects. Bulk tissue cannot separate tumor-intrinsic change from cell-composition change; external single-cell and spatial datasets used for context were not from these patients. The GSE155434 transcript comparison is small (11 pairs) and cross-cohort. Survival-trajectory associations carry immortal-time/landmark bias. No perturbation experiments were performed.

## 4.6 Future directions

Priorities: (i) same-patient longitudinal RNA integration (PRJNA1051047; 70 patients with complete triple-omics pairing are available for processed-matrix analysis) to directly quantify transcript–protein coupling over recurrence; (ii) spatial deconvolution of the inflammatory–metabolic axis (GeoMx cohorts) to determine its cellular source; (iii) perturbation tests of the class-B metabolic/inflammatory programs in GBM models across proliferation states.

# 5. Conclusions

GBM recurrence follows a two-axis, patient-resolved evolutionary structure: a tightly coupled proliferative-machinery contraction and a partially independent inflammatory–metabolic up-shift. Most longitudinal proteome and phosphoproteome changes are state-coupled or abundance-driven. Competing-hypothesis adjudication — not any single molecule — is the central deliverable: module-specific recurrence claims should be required to survive it.

# Acknowledgements

The authors thank the data-generating consortia and study authors (CPTAC, TCGA, GLASS, DepMap, GDSC, GEO contributors, and the KNCC/ICPC investigators of the reference longitudinal study).

# Authors' contributions

JC: conceptualization, data curation, formal analysis, methodology, software, validation, visualization, writing – original draft. JY, BY, TZ, PZ: validation, writing – review & editing.

# Ethics approval and consent to participate

Only publicly available, de-identified datasets were analyzed under each provider's terms. No new human data were collected.

# Consent for publication

Not applicable.

# Data availability

PDC000514/PDC000515 (Proteomic Data Commons); CPTAC-GBM (LinkedOmics); GSE155434, GSE254875 (GEO); TCGA (UCSC Xena/GDC); TCGASpliceSeq (MD Anderson); DepMap 24Q2; GDSC2; OmniPath; ComplexPortal. Code: https://github.com/huozheshuashuai/GBM_WEBS. The patient-level data dictionary accompanies the code release.

# Competing interests

The authors declare no competing interests.

# Funding

No specific funding was received for this work.

# References

[to be completed — key: Kim KH et al. Cancer Cell 2024; CPTAC-GBM confirmatory study; TCGASpliceSeq; DepMap; GDSC2; OmniPath; ComplexPortal; GSE155434 study; Neftel et al. 2019]

# Figure Legends

**Figure 1. Study design and error-aware QC.** (a) Cohort structure: 105 paired proteomes, 89 paired phosphoproteomes, 70 triple-omics patients; 100/105 pairs within the same TMT plex (batch-controlled deltas by design). (b) Technical variance from cross-plex duplicate aliquots vs inter-patient module-delta variance (20–260× ratios).

**Figure 2. Patient-resolved evolution map.** (a) Module-delta co-evolution heatmap (23 modules, per-patient Spearman). (b) Two-axis structure (PC1 41%, PC2 38%); patients colored by proliferation delta. (c) Bootstrap CIs of top co-evolving pairs. (d) Silhouette analysis: continuum, no stable discrete subtypes.

**Figure 3. Competing-hypothesis adjudication.** (a) Module screen: raw vs state-adjusted effects; class A (state-explained) vs class B (independent residual, all up, metabolic/inflammatory). (b) Adjudication model schematic and coefficient summary. (c) Proliferation-matched tertiles for the spliceosome.

**Figure 4. Worked example: the spliceosome and matched complex controls.** (a) Paired volcano (105/110 down). (b) Subcomplex uniformity (KW p = 0.29). (c) Adjusted intercepts across three proliferation definitions. (d) Matched complex controls: machinery contraction vs spared proteasome. (e) No transcript-level counterpart (GSE155434).

**Figure 5. Phosphoproteome: abundance shadows.** (a) Coverage audit (pairs/site distribution). (b) Raw vs parent-corrected site deltas. (c) KSEA-style kinase analysis: 5 raw hits, 0 after correction, LOO sensitivity. (d) Module quadrant: protein-only changes dominate.

**Figure 6. External replication, prediction, clinical association, and boundaries.** (a) CPTAC cross-sectional replication of the KNCC module co-variation structure (21 module pairs, ρ = 0.94, p = 5 × 10⁻⁶). (b) M0/M1/M2 cross-validated error with permutation p-values (no gain from broader baseline state). (c) OS by proliferation-trajectory tertile (p = 0.004; exploratory, landmark caveats). (d) Boundary summary: DepMap, GDSC2, CPTAC corroboration of PTM stoichiometry.

# Tables

**Table 1.** Cohorts, data layers, and patient-level intersection.

**Table 2.** Module-level adjudication (50 Hallmark modules: raw effect, adjusted intercept, class).

**Supplementary Tables.** Per-protein paired deltas; per-site phospho results (raw and corrected); kinase-substrate results with LOO; complex-control matched-null results; prediction model details.
