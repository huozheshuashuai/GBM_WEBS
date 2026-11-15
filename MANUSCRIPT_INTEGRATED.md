# Title

Protein-level regulation of the spliceosome across the glioblastoma state spectrum and recurrence: an integrated cross-sectional and longitudinal proteogenomic analysis

# Running title

Spliceosome proteome dynamics in GBM recurrence

# Authors

Jia Chen*, Junhan Yang, Bing Yang, Tianyi Zhou, Peijun Zhao

Department of Bioinformatics, Chengdu Neusoft University, Chengdu, China

\* Corresponding author: Jia Chen. Email: 1176625889@qq.com. ORCID: 0009-0002-0082-3050

# Abstract

**Background.** Spliceosomal dysregulation is a recognized feature of cancer, yet how the spliceosome itself is regulated across glioblastoma (GBM) tumors and over recurrence remains poorly characterized. Transcript- and protein-level studies have yielded contradictory impressions, and no longitudinal proteomic map of the spliceosome in GBM exists.

**Methods.** We integrated single-cell and bulk transcriptomics with the CPTAC-GBM proteome, acetylome and phosphoproteome (n = 99), the longitudinal KNCC proteome and phosphoproteome (105 paired primary–recurrence proteomes; 89 patients with paired phosphoproteomes), TCGA SpliceSeq alternative-splicing profiles (TCGA-GBM and TCGA-LGG), and paired recurrent RNA-seq (GSE155434). Per-protein mRNA-adjusted residuals quantified protein–transcript discordance; paired models with competing-hypothesis adjustment separated cell-state-coupled from independent changes.

**Results.** Across CPTAC tumors, inter-patient spread of spliceosome protein–mRNA discordance exceeded abundance- and variance-matched random complexes (10,000 permutations, empirical p < 10⁻⁴). Longitudinally, 105 of 110 spliceosome proteins decreased at recurrence (paired Wilcoxon; direction-uniform), with no subcomplex selectivity (Kruskal–Wallis p = 0.29). Formal adjudication attributed 80–85% of the decline to the proliferative-to-neuronal state transition, leaving a small residual (−0.04 log2) that was consistent across three proliferation definitions, patient bootstrap, and leave-one-patient-out analysis; patients whose tumors gained proliferation at recurrence showed the opposite (spliceosome increase, p = 0.001). Recurrence-associated spliceosome protein decline had no transcript-level counterpart in independent paired RNA-seq (p = 0.77), indicating protein-level regulation. Phosphorylation and acetylation changes tracked protein abundance in both cohorts; after parent-protein correction, modification density was reduced, and no abundance-independent spliceosome phosphosite change survived. Spliceosome-machinery expression tracked 1,484 alternative splicing events in TCGA-GBM, replicating at 72% sign concordance in TCGA-LGG (ρ = 0.44). Proliferation-trajectory subgroups differed in overall survival (log-rank p = 0.004; exploratory, subject to landmark considerations). Dependency and drug-sensitivity screens on the spliceosome axis were negative at FDR thresholds and are reported as boundary results.

**Conclusions.** The GBM spliceosome is remodelled predominantly at the protein level, both across the inter-tumor state spectrum and during recurrence. Most longitudinal change is coupled to the cell-state transition, with a small robust independent component and altered splicing output. These data define the boundaries of spliceosome regulation in GBM and argue against spliceosome-intrinsic reprogramming as a primary recurrence driver.

**Keywords:** glioblastoma; spliceosome; proteogenomics; longitudinal proteomics; alternative splicing; recurrence; protein–mRNA discordance

# 1. Introduction

Glioblastoma (GBM) recurs nearly universally despite maximal therapy, and recurrent tumors differ from their primary counterparts in proliferation, neuronal-like programs, and microenvironment composition. The spliceosome is essential for transcriptome maturation and is increasingly implicated in cancer, but reports disagree on whether spliceosomal activity in glioma is elevated, suppressed, or prognostically favorable. A recurring difficulty is that transcript-level measurements do not necessarily reflect protein abundance of the splicing machinery, and post-translational regulation of spliceosome proteins is well documented in other systems.

Two designs are needed to resolve this. Cross-sectional proteogenomic cohorts can test whether spliceosome protein levels track their transcripts across patients. Longitudinal paired cohorts can test whether the spliceosome is remodelled between primary and recurrent disease within the same patient. Neither design has been applied to the spliceosome in GBM at proteomic depth. The longitudinal CPTAC/KNCC cohort (Kim et al., 2024, Cancer Cell) provides 105 paired primary–recurrence proteomes, and its companion study of alternative-splicing evolution (GSE155434) provides paired RNA; the CPTAC-GBM confirmatory cohort provides paired RNA–protein–PTM data for 99 tumors.

We asked four questions: (i) Is spliceosome protein abundance coupled to its mRNA across GBM tumors? (ii) Is the spliceosome remodelled at recurrence, and is any remodeling selective to particular subcomplexes? (iii) How much of any change is explained by the proliferative-to-neuronal state transition, and does an independent component remain? (iv) Does spliceosome variation have measurable consequences for alternative splicing and clinical outcome? We report positive and negative results explicitly, including boundary analyses that constrain mechanistic interpretation.

# 2. Methods

## 2.1 Data sources

Cross-sectional: CPTAC-GBM confirmatory cohort (LinkedOmics; TMT proteome 11,141 proteins, acetylome, phosphoproteome, RNA-seq; n = 99 patients with paired RNA–protein). Longitudinal: KNCC GBM-evolution study (PDC000514 proteome, PDC000515 phosphoproteome; clinical from the published Table S1). Splicing: TCGASpliceSeq PSI for TCGA-GBM (28,071 events, 153 samples) and TCGA-LGG (515 samples). Paired recurrent RNA-seq: GSE155434 (37 samples; exon-level counts; 11 patients with primary–recurrent pairs). Cell models: DepMap 24Q2 (CRISPR gene effect, expression; 126 brain-lineage lines with both) and GDSC2 (fitted dose responses). Single-cell and spatial datasets from the prior study (GSE138794, GSE131928, GSE194329) were used for state-annotation support only.

## 2.2 Proteome processing and quality control (KNCC)

PDC processed TMT11 log2 ratios (vs pooled reference) were used. "Other"/reference channels were excluded; duplicate aliquots (same sample in two plexes) were averaged; 11,320 proteins × 214 tumor samples (109 patients; 105 with paired T1/T2) were retained. Median per-sample missingness was 7.4%. PCA was used to inspect structure; PC1 associated with timepoint (p = 1.5 × 10⁻⁹) and PC2 with TMT plex (p = 2.8 × 10⁻⁴), motivating paired designs with plex-aware sensitivity analyses.

## 2.3 Discordance axis (CPTAC)

Per-protein residuals from per-patient linear models (protein ~ mRNA) quantify protein abundance beyond transcript prediction. The spliceosome protein discordance score (SPDS) is the per-patient mean residual across the 109 quantified spliceosome proteins. Specificity was tested against 10,000 random gene sets matched for abundance decile and within-decile variance tercile, drawn without replacement from non-spliceosome genes (empirical p).

## 2.4 Longitudinal and competing-hypothesis models (KNCC)

Per-patient paired deltas (T2 − T1) were computed per protein. Per-protein paired Wilcoxon tests with Benjamini–Hochberg correction assessed recurrence change. The competing-hypothesis model regressed per-patient module deltas on proliferation (14-marker core; sensitivity: Hallmark E2F targets and G2M checkpoint), neuronal-marker and immune-marker deltas (uncentered predictors; the intercept is the expected module change when state scores do not change). VIF, partial R², patient bootstrap (2,000×) and leave-one-patient-out analysis assessed robustness. Proliferation-matched tertile analysis tested whether spliceosome change persists within strata of proliferation change. Subcomplex assignments used ComplexPortal (9606) complex names with overlapping membership allowed; module selectivity was tested by Kruskal–Wallis across per-protein deltas.

## 2.5 Phosphoproteome and acetylome

KNCC phosphosite deltas (89 paired patients) were tested per site (paired Wilcoxon; ≥25 pairs required). Parent-protein correction used two methods: delta subtraction and per-site regression of site delta on parent-protein delta (residuals; primary). CPTAC site-level acetylation/phosphorylation were correlated with SPDS (Spearman), raw and after parent-protein correction. Modification-density interpretation used corrected values only.

## 2.6 Alternative splicing

Per-event Spearman correlations between PSI and the spliceosome mRNA score (mean expression of KEGG spliceosome genes) were computed in TCGA-GBM and TCGA-LGG; events with q < 0.05 and |ρ| > 0.2 in GBM were tested for sign concordance and effect-size correlation in LGG. GSE155434 exon counts were summed to gene level for state-score analyses.

## 2.7 Clinical analyses

KNCC clinical variables (PFS, OS, vital status) were linked by PDC_ID. Cox models adjusted for age and sex. Proliferation-trajectory groups (tertiles of proliferation-score delta) were compared by log-rank tests. Because trajectory assignment requires survival to recurrence, these analyses are exploratory and subject to immortal-time/landmark considerations.

## 2.8 Dependency and pharmacogenomic screens

Spliceosome expression scores in DepMap brain-lineage lines were correlated with CRISPR gene effects of spliceosome genes (Spearman, BH-FDR) and with GDSC2 fitted dose responses (LN_IC50) across 271 drugs.

## 2.9 Statistics and reproducibility

Analyses used R 4.5.2. Paired designs were used throughout the longitudinal analyses; patients, not samples, were the statistical unit. Permutation controls excluded the tested gene set from the background pool, matched abundance and variance, and sampled without replacement. The full pipeline (scripts 86–100) is available at https://github.com/huozheshuashuai/GBM_WEBS. During preparation, generative AI tools were used for language editing and code-formatting assistance; all analyses were executed and verified by the authors.

# 3. Results

## 3.1 Inter-patient spliceosome protein–mRNA discordance is specific and tracks the tumor state spectrum

In 99 CPTAC-GBM tumors with paired RNA and protein, spliceosome complex-level protein and mRNA scores were correlated (ρ = 0.66), so the layers are not globally decoupled. However, mRNA explained only part of protein variance: per-protein mRNA-adjusted residuals showed systematic positive loading on a per-patient discordance axis (SPDS), and the inter-patient spread of SPDS exceeded 10,000 abundance- and variance-matched random gene sets (empirical p < 10⁻⁴; Fig. 1). Non-spliceosome proteins tracking SPDS were dominated by neuronal markers in the opposite direction (NEFL, NEFH, NEFM, STX1A/1B, STXBP1, CAMK2A, SNAP25 family), indicating that the discordance axis is anchored in the proliferation-versus-neuronal-content state spectrum rather than being spliceosome-autonomous.

## 3.2 The spliceosome proteome declines at recurrence without subcomplex selectivity

In 105 paired KNCC primary–recurrence proteomes, 105 of 110 quantified spliceosome proteins decreased at recurrence (FDR < 0.05; direction-uniform: 109/110 negative), with a median paired delta of −0.20 log2 against a genome-wide median of −0.006 (Wilcoxon p = 2.8 × 10⁻⁴⁵); the complex-level score declined in paired analysis (p = 3.5 × 10⁻⁸), exceeding matched random sets (10,000 strict permutations, empirical p = 10⁻⁴; Fig. 2). However, decline was uniform across U1, U2, U4/U6, U5, tri-snRNP, NineTeen/Prp19, SR-protein, hnRNP and DDX/DHX modules (ComplexPortal-annotated; Kruskal–Wallis p = 0.29): there is no evidence for selective subcomplex remodeling.

## 3.3 Most of the decline is coupled to the cell-state transition; a small independent component remains

Formal competing-hypothesis models regressed per-patient spliceosome deltas on proliferation, neuronal and immune marker deltas (no gene overlap between marker sets and the spliceosome; VIF < 2.2). Across three independent proliferation definitions, 80–85% of the decline was attributable to the state transition, leaving a small adjusted intercept of −0.03 to −0.04 log2 (bootstrap 95% CI excluding zero; leave-one-patient-out stable; Fig. 3). Proliferation-stratified analysis showed the coupling directly: patients with falling proliferation showed the steepest spliceosome decline (−0.42), whereas patients whose tumors gained proliferation at recurrence showed a spliceosome increase (+0.10, p = 0.001). Independent paired RNA-seq (GSE155434, 11 pairs) showed no transcript-level spliceosome decline at recurrence (median +0.03, p = 0.77), consistent with protein-level rather than transcriptional regulation, although the small sample limits power. In a whole-proteome extension (50 Hallmark modules), 20 modules were fully explained by the state transition while 11 retained small independent residuals — all positive (up at recurrence) and dominated by metabolic and inflammatory programs (reactive oxygen species pathway, IL2–STAT5 signaling, xenobiotic/fatty-acid metabolism).

## 3.4 Post-translational modification changes are stoichiometric, not independent

Spliceosome phosphosites declined at recurrence in raw KNCC data (47/107 testable sites FDR < 0.05, all down), but after regression-based parent-protein correction none remained significant (0/107), and the same pattern held genome-wide at module level. In CPTAC, raw acetylation and phosphorylation of spliceosome proteins tracked the discordance axis (OR ≈ 9.7 for both), but after parent-protein correction the direction inverted (acetylation OR = 0.43; phosphorylation OR = 0.19; both significant), i.e., per-copy modification density is lower on the more abundant spliceosome proteins. The two cohorts therefore agree that spliceosome PTM changes are largely consequences of abundance change (stoichiometric dilution), not independent regulatory remodeling (Fig. 4).

## 3.5 Spliceosome variation is accompanied by altered splicing output

In TCGA-GBM, 1,484 of 28,071 alternative-splicing events tracked the spliceosome mRNA score (q < 0.05, |ρ| > 0.2). Of the 1,425 testable in TCGA-LGG, 72.2% were sign-concordant (effect-size ρ = 0.44, p = 3 × 10⁻⁶⁹; Fig. 5). Top events included known splicing-regulatory genes (e.g., an RBM4 exon-skipping event), consistent with feedback within the splicing network. This replication is weaker than transcriptome-wide axes but supports a measurable splicing-output correlate of machinery abundance.

## 3.6 Clinical trajectory heterogeneity and boundary results

Proliferation-trajectory subgroups (drop/stable/rise tertiles of the proliferation-score delta) differed in overall survival (median OS 870/730/945 days; log-rank p = 0.004) and progression-free survival (p = 0.044; Fig. 6). Because trajectory assignment requires survival to recurrence, these associations are exploratory and subject to immortal-time considerations. Spliceosome scores themselves (T1 level or longitudinal delta) were not prognostic after age/sex adjustment (OS p = 0.31; PFS p = 0.28). Dependency and drug screens on the spliceosome axis were negative: no spliceosome gene showed FDR-significant score-dependent essentiality in DepMap brain lines, and no GDSC2 compound met FDR < 0.1. We report these as boundary results that constrain therapeutic extrapolation.

# 4. Discussion

## 4.1 Principal findings

This study provides, to our knowledge, the first paired longitudinal proteomic map of the spliceosome in GBM and a cross-sectional discordance analysis in an independent cohort. Four conclusions emerge. First, spliceosome protein abundance varies across patients partly independently of its transcripts, and the spread of this discordance exceeds matched controls. Second, the spliceosome proteome declines at recurrence, uniformly across subcomplexes. Third, the large majority of this decline is explained by the proliferative-to-neuronal state transition — the same axis that anchors the cross-sectional discordance — with a small robust residual and a protein-level (not transcriptional) longitudinal signature. Fourth, spliceosome PTM changes in both cohorts are largely stoichiometric consequences of abundance, and machinery variation is accompanied by modest, replicable changes in splicing output.

## 4.2 Relationship to prior work

Kim et al. (2024) reported the proliferative-to-neuronal transition in this longitudinal cohort; our analysis quantifies how much of the spliceosome change that transition explains and adds the competing-hypothesis adjudication, the subcomplex-selectivity test, the PTM stoichiometry analysis, and the splicing-output correlation. The paired-splicing study (GSE155434) reported alternative-splicing evolution under therapy; we add the negative but informative observation that the spliceosome machinery itself does not decline transcriptionally in that cohort. Our results do not support spliceosome-intrinsic reprogramming as a primary driver of recurrence; rather, spliceosome abundance is a tightly coupled readout of the state transition, with possible small independent components whose biological significance remains to be established.

## 4.3 Methodological contribution

The paired competing-hypothesis framework used here — uncentered state-covariate models with intercept interpretation, marker-overlap audit, multiple independent state scores, patient bootstrap, leave-one-patient-out analysis, matched-strata analysis, and strict matched permutation controls that exclude the tested set — offers a template for proteogenomic claims of module-specific regulation. We recommend that claims of selective complex remodeling be required to survive this class of adjudication, since prominent candidate modules (here, the spliceosome) can otherwise be credited with changes that are generic consequences of cell-state transitions.

## 4.4 Limitations

KNCC RNA-seq (PRJNA1051047) was not analyzed here; same-patient protein–RNA longitudinal pairing would strengthen the protein-level-regulation claim beyond the cross-cohort consistency we present. The GSE155434 mRNA analysis is underpowered (11 pairs) and cannot exclude modest transcript-level change. Phosphoproteome missingness (median 74% of sites per patient-pair) limits site-level power, and TMT ratio compression may attenuate true fold changes. Clinical trajectory analyses are subject to immortal-time/landmark bias. Proliferation/neuronal markers are protein-level proxies, not direct cell-composition measurements; single-cell or spatial deconvolution of the longitudinal cohort would refine attribution. Finally, all analyses are observational; no perturbation experiments were performed, and no therapeutic implications are claimed.

## 4.5 Future directions

Priority experiments include: (i) same-patient longitudinal RNA-seq integration (PRJNA1051047) to directly test transcript–protein coupling over recurrence; (ii) targeted proteomics of candidate spliceosome proteins across the state spectrum; (iii) perturbation of spliceosome components in GBM models across proliferation states to test whether the small independent residual is causal; and (iv) deeper investigation of the independent metabolic/inflammatory programs (ROS, IL2–STAT5) that emerged from the whole-proteome screen.

# 5. Conclusions

The GBM spliceosome is remodelled predominantly at the protein level across the inter-tumor state spectrum and during recurrence. Most longitudinal change is coupled to the proliferative-to-neuronal transition; a small independent component and altered splicing output remain. Rigorous competing-hypothesis adjudication, rather than any single molecular claim, is the central deliverable of this study.

# Acknowledgements

The authors thank the data-generating consortia (CPTAC, TCGA, GLASS, DepMap, GDSC, GEO contributors, and the KNCC/ICPC investigators) and the original study authors.

# Authors' contributions

JC: conceptualization, data curation, formal analysis, methodology, software, validation, visualization, writing – original draft. JY, BY, TZ, PZ: validation, writing – review & editing.

# Ethics approval and consent to participate

Only publicly available, de-identified datasets were analyzed, in accordance with each provider's data-use terms. No new human data were collected.

# Consent for publication

Not applicable.

# Data availability

All datasets are public: CPTAC-GBM (LinkedOmics), PDC000514/PDC000515 (Proteomic Data Commons), GSE155434 and GSE254875 (GEO), TCGASpliceSeq (MD Anderson), TCGA cohorts (UCSC Xena/GDC), DepMap 24Q2, GDSC2. Analysis code: https://github.com/huozheshuashuai/GBM_WEBS.

# Competing interests

The authors declare no competing interests.

# Funding

No specific funding was received for this work.

# References

[to be completed on revision — key citations: Kim KH et al. Cancer Cell 2024 (PDC000514/515); the GSE155434 paired-splicing study; CPTAC-GBM confirmatory proteogenomics; TCGASpliceSeq (Ryan et al.); DepMap (Tsherniak et al.); GDSC2 (Yang et al.); ComplexPortal; KEGG spliceosome; Neftel et al. Cell 2019; clinical and methods references as in the prior version]

# Figure Legends

**Figure 1. Cross-sectional spliceosome protein–mRNA discordance (CPTAC-GBM, n = 99).** (a) Complex-level spliceosome protein score vs mRNA score (ρ = 0.66). (b) Inter-patient variance of the discordance score vs 10,000 abundance- and variance-matched random gene sets (empirical p < 10⁻⁴). (c) Non-spliceosome proteins tracking the discordance axis are dominated by neuronal markers (opposite direction).

**Figure 2. Longitudinal spliceosome decline (KNCC, 105 pairs).** (a) Paired per-protein deltas (T2 − T1) with spliceosome highlighted (105/110 FDR < 0.05, all down). (b) Complex score by timepoint (paired Wilcoxon p = 3.5 × 10⁻⁸). (c) Strict matched-permutation null (empirical p = 10⁻⁴). (d) Subcomplex uniformity (ComplexPortal modules; Kruskal–Wallis p = 0.29).

**Figure 3. Competing-hypothesis adjudication.** (a) Per-patient spliceosome delta vs proliferation delta (ρ = 0.79). (b) Adjusted intercepts across three proliferation definitions with bootstrap CIs. (c) Proliferation-matched tertiles: decline in falling-proliferation patients, increase in rising-proliferation patients. (d) Whole-proteome module screen: 20 modules state-explained, 11 with small independent residuals (all up-regulated; metabolic/inflammatory).

**Figure 4. PTM stoichiometry.** (a) KNCC raw vs parent-corrected spliceosome phosphosite deltas. (b) CPTAC spliceosome acetylation/phosphorylation vs discordance axis, raw vs corrected (direction inversion after correction). (c) Modification-density dilution schematic.

**Figure 5. Splicing output.** (a) PSI vs spliceosome mRNA score in TCGA-GBM (1,484 events). (b) Cross-cohort replication in TCGA-LGG (72% sign-concordant; ρ = 0.44). (c) Event-type composition and top events (RBM4, ILF3, PNKD, NDE1).

**Figure 6. Clinical trajectories.** (a) Overall survival by proliferation-trajectory group (log-rank p = 0.004; exploratory, landmark caveats). (b) PFS (p = 0.044). (c) Null results summary: spliceosome scores vs survival; DepMap dependency and GDSC2 sensitivity screens (no FDR-significant hits).

# Tables

**Table 1. Cohorts and data layers.**

**Table 2. Module-level adjudication results (50 Hallmark modules): class, raw and adjusted effects.**

**Supplementary Tables.** Per-protein paired deltas (KNCC); per-site PTM results (raw and corrected, both cohorts); per-event PSI correlations (GBM, LGG); GDSC2 and DepMap full result tables.
