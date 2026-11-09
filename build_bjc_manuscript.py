import re

s = open("MANUSCRIPT_GND.md", encoding="utf-8").read()

# ---------- 1. neutral descriptive title (BJC requirement) ----------
s = s.replace("**Title**: A writer-dominant expression state is associated with post-transcriptional spliceosome accumulation in glioblastoma",
              "**Title**: Writer–eraser balance and post-transcriptional spliceosome accumulation in glioblastoma")

# ---------- 2. 3.5 section title softening ----------
s = s.replace("## 3.5 Mechanistic evidence: splicing-factor acetylation couples to protein abundance",
              "## 3.5 Post-translational modification patterns associated with spliceosome accumulation")

# ---------- 3. compression: Introduction (~130 words) ----------
old_intro2 = """What upstream program could drive such layer-specific accumulation? Histone acylation — including acetylation and lactylation of histone tails — is written by histone acetyltransferases (EP300, KAT2A, KAT2B, KAT5, KAT8, CREBBP) and erased by class I HDACs and sirtuins (HDAC1–3, SIRT1–3, SIRT6–7) [5,10–12], and the net modification state reflects the balance between these opposing enzymatic capacities. Acylation is disease-relevant in glioma [6–9,13–15], and several lactylation-associated prognostic signatures have been derived [27–32], but these remain at the transcriptomic correlation level. Importantly, the same enzymes also catalyze non-histone protein acetylation, including of splicing factors — a regulatory layer distinct from chromatin-level histone modification. We therefore hypothesized that writer–eraser imbalance could stabilize the spliceosome at the protein level without altering its transcription — a mechanism that may not be fully captured by transcriptomic measurements."""
new_intro2 = """What upstream program could drive such accumulation? Histone acylation is written by acetyltransferases (EP300, KAT2A, KAT2B, KAT5, KAT8, CREBBP) and erased by class I HDACs and sirtuins (HDAC1–3, SIRT1–3, SIRT6–7) [5,10–12], and these enzymes also acetylate non-histone proteins, including splicing factors. We therefore hypothesized that writer–eraser imbalance could stabilize the spliceosome at the protein level without altering its transcription — a mechanism that may not be fully captured by transcriptomic measurements."""
s = s.replace(old_intro2, new_intro2)

# ---------- 4. compression: Methods ----------
# 2.1 data sources (compress enumeration)
s = s.replace(
"scRNA-seq: GSE138794 (35,823 QC-passed cells, 20 samples, 17 patients). Spatial: GSE194329 (six GBM Visium slices). External scRNA: GSE131928 (Smart-seq2 and 10x). Bulk: TCGA-GBM (LinkedOmics legacy, n = 144; GDC STAR-counts, n = 371), CGGA-693, CGGA-325, GSE7696, GSE83300 (n = 50); RPPA: TCPA RPPA500 (n = 71 matched). Longitudinal: GLASS consortium RNA-seq (gene-level TPM, 425 aliquots; 403 mapped to cases; 107 primary tumors with survival; 179 primary–recurrence pairs from 148 patients; 213 recurrence-interval records from 173 patients) accessed via Synapse (syn17038081) [38]. Proteomics: CPTAC-GBM TMT proteome [24] (11,141 proteins), acetylome (13,556 sites), phosphoproteome (31,010 quantified sites), and RNA-seq (n = 99). Alternative splicing: TCGASpliceSeq percent-spliced-in (PSI) values for TCGA-GBM (28,071 events, 153 samples matched to WEBS) [39]. Pharmacogenomics: DepMap 24Q2 [25] (CRISPRGeneEffect, expression, model metadata) and GDSC2 [26] (fitted dose responses). Pan-cancer: all 33 TCGA cohorts (STAR-TPM and survival, UCSC Xena GDC hub). All cohorts are summarized in Table 1.",
"Data sources (all summarized in Table 1): scRNA-seq GSE138794 (35,823 cells, 17 patients); spatial GSE194329; external scRNA GSE131928; bulk TCGA-GBM (LinkedOmics n = 144; GDC n = 371), CGGA-693/325, GSE7696, GSE83300; longitudinal GLASS (425 aliquots; 107 primary with survival; 179 pairs/148 patients; 213 recurrence-interval records/173 patients) [38]; CPTAC-GBM proteome/acetylome/phosphoproteome/RNA-seq (n = 99) [24]; TCGASpliceSeq PSI for TCGA-GBM (28,071 events, 153 samples) [39]; DepMap 24Q2 [25]; GDSC2 [26]; all 33 TCGA cohorts (UCSC Xena GDC hub).")

# 2.3 single-cell (tighten)
s = s.replace(
"Seurat v5 pipeline [16] (log-normalization, 2,000 variable features, PCA dims 1:30, UMAP, resolution 0.8). Cell identities were anchored to the original study's annotations and cross-checked by marker-based scoring. CopyKAT (per-sample, ≤1,200 cells, hg19) distinguished aneuploid tumor cells. Region types: Core (MKI67, TOP2A, PCNA), Margin (FN1, ITGA5, SPP1, MMP2), Normal (RBFOX3, MAP2, GFAP) composite max-score assignment. Hypoxia decoupling: logistic regression (Margin ~ WEBS + Hallmark Hypoxia) and partial correlation. Trajectory: Slingshot [18] on neoplastic cells with LLS origin; GAM fit of WEBS along pseudotime. Spliced/unspliced matrices are not deposited for GSE138794, precluding RNA velocity. Pseudoreplication control: key single-cell associations were re-tested with patient as the statistical unit — patient × state pseudobulk means, linear mixed-effects models with patient as random intercept (lme4/lmerTest), and per-patient Spearman correlations with one-sample tests (Fig. S21).",
"Seurat v5 pipeline [16]; identities anchored to original annotations; CopyKAT distinguished aneuploid tumor cells. Region types (Core/Margin/Normal-like) by composite marker scores. Hypoxia decoupling by logistic regression and partial correlation. Trajectory by Slingshot [18] (LLS origin; RNA velocity unavailable). Pseudoreplication control: key associations re-tested with patient as the statistical unit (patient × state pseudobulk means; mixed-effects models with patient random intercept; per-patient correlations; Fig. S21).")

# 2.4 + 2.5 (tighten)
s = s.replace(
"SCTransform per slice; RCTD was attempted [19] but its internal caching step failed on our platform, so label-transfer-free spot-level WEBS was used; spot-level WEBS on SCT-normalized data; dual-hotspots defined at the 75th percentile of writer and metabolism scores; Moran's I and Geary's c with six-nearest-neighbor weights and 1,000 label permutations (spdep).",
"SCTransform per slice; spot-level WEBS; Moran's I and Geary's c with 1,000 label permutations (RCTD failed on our platform [19]).")

# 2.6 proteomics (compress filter prose)
s = s.replace(
"limma differential abundance (LHS vs LLS) on the CPTAC TMT proteome; spliceosome enrichment by Fisher's exact test (111 spliceosome proteins passing detection filtering out of the curated list). Layer discordance: protein WEBS vs protein pathway score, and mRNA WEBS vs mRNA pathway score (Spearman), with patient-paired bootstrap 95% CIs (2,000 resamples of patients, RNA and protein drawn jointly per patient), leave-one-protein-out analysis, partial correlations controlling proliferation/immune/mesenchymal/hypoxia/OPC protein programs, and specificity against 10,000 random gene sets matched for proteome abundance (empirical p). Acetylome: 13,556 measured sites were filtered to 8,350 sites quantified in ≥50% of samples for limma testing; enrichment of spliceosome sites among hyper-acetylated sites; per-protein acetylation–abundance coupling; site-level ordered trend across LLS/LDTS/LHS (Jonckheere trend test). Phosphoproteome: per-site limma (LHS vs LLS; sites quantified in ≥50% of samples); Fisher enrichment of spliceosome-protein phosphosites among significant sites; SRSF-family phosphosites examined separately; kinase abundance (SRPK1/2, CLK1–4, CDK7/9/12/13) correlated with WEBS (Spearman). To separate occupancy-level from abundance-level changes, PTM analyses were repeated on abundance-corrected site values (site intensity − parent-protein abundance, log scale). Histone acetylation: mean across histone H1/H2A/H2B/H3/H4 sites.",
"limma differential abundance (LHS vs LLS) on the CPTAC proteome; spliceosome enrichment by Fisher's test (111 proteins). Layer discordance (protein vs mRNA pathway scores vs respective WEBS, Spearman) with patient-paired bootstrap CIs (2,000 resamples), leave-one-protein-out, partial correlations controlling proliferation/immune/MES/hypoxia/OPC programs, and 10,000 expression-matched random gene sets (empirical p). Acetylome (13,556→8,350 sites quantified in ≥50% of samples) and phosphoproteome (31,010 sites): per-site limma; spliceosome enrichment; per-protein coupling; SRSF-family and kinase analyses (SRPK1/2, CLK1–4, CDK7/9/12/13). PTM analyses were repeated on abundance-corrected site values (site − parent protein, log scale). Histone acetylation: mean across histone sites.")

# ---------- 5. compression: Results (remove duplicated numbers) ----------
# 3.1 (tighten)
s = s.replace(
"We assembled 35,823 cells passing quality control from 20 tumor samples of 17 GBM patients (GSE138794); CopyKAT-purified malignant cells (used for trajectory, region, and pseudobulk analyses) derived from 8 of these patients. WEBS varied continuously across cells and differed between Neftel malignant states in both the discovery cohort and an external cohort (GSE131928; Smart-seq2: Kruskal–Wallis p = 1.6 × 10⁻¹⁶⁸; 10x: p = 1.0 × 10⁻³¹³), with OPC-like cells lowest (Fig. 2, Fig. S25; external-cohort panels in Fig. S24). Copy-number-based tumor purification (CopyKAT; 16,255 cells called) confirmed that stratification is preserved within 4,421 aneuploid cells at the cell level (p < 10⁻¹⁶), excluding normal-cell admixture as an explanation (Fig. S1); note that this cell-level statistic does not by itself establish patient-level significance, which is addressed separately below. To guard against pseudoreplication (cells from the same patient are not independent replicates), all key single-cell associations were re-tested with patient as the statistical unit: the Neftel-state gradient remained significant under pseudobulk mixed-effects modeling (p = 7.4 × 10⁻⁴), and the within-malignant-cell hypoxia decoupling was confirmed per-patient (median ρ = −0.02, one-sample Wilcoxon p = 1), whereas the regional gradient and the pseudotime trend attenuated at patient level (p = 0.77 and p = 0.22, respectively) and are therefore presented as supportive cell-level observations (Fig. S21).",
"We assembled 35,823 cells from 17 GBM patients (GSE138794). WEBS varied continuously and differed between Neftel malignant states in both the discovery cohort and an external cohort (GSE131928), with OPC-like cells lowest (Fig. 2, Figs. S24–S25). CopyKAT purification confirmed stratification within 4,421 aneuploid cells (p < 10⁻¹⁶; Fig. S1). To guard against pseudoreplication, all key associations were re-tested with patient as the statistical unit: the Neftel-state gradient remained significant (pseudobulk mixed-effects p = 7.4 × 10⁻⁴) and hypoxia decoupling was confirmed per-patient, whereas the regional gradient and pseudotime trend attenuated at patient level and are presented as supportive observations (Fig. S21).")

# 3.2 (tighten)
s = s.replace(
"WEBS differed across Core-like, Margin-like, and Normal-like regions when all cells were considered (Kruskal–Wallis p < 2.2 × 10⁻¹⁶; Fig. 4a), and a patient-level Margin-vs-Core paired test on all cells confirmed the regional gradient (paired Wilcoxon, n = 17; Fig. 4b). However, within CopyKAT-purified malignant cells the regional gradient attenuated at the patient level (pseudobulk mixed-effects p = 0.77; only 3 patients contributed both Margin-like and Core-like malignant cells), so the regional association is reported as a whole-tumor-compartment observation rather than a malignant-cell-intrinsic one. Although WEBS correlated with the Hallmark Hypoxia signature across all cells, this association was driven by non-malignant compartments: within malignant cells WEBS was essentially decoupled from hypoxia (ρ = −0.054; confirmed per-patient, median ρ = −0.02, one-sample Wilcoxon p = 1), and the Margin association persisted after hypoxia adjustment (logistic regression p = 3.8 × 10⁻⁴; Fig. 4c,d).",
"WEBS differed across Core-like, Margin-like, and Normal-like regions (all cells, p < 2.2 × 10⁻¹⁶; patient-paired n = 17; Fig. 4a,b), although within purified malignant cells the gradient attenuated at the patient level (p = 0.77). WEBS was decoupled from hypoxia within malignant cells (ρ = −0.054; per-patient median ρ = −0.02, p = 1), and the Margin association persisted after hypoxia adjustment (p = 3.8 × 10⁻⁴; Fig. 4c,d).")

# 3.3 (tighten)
s = s.replace(
"Slingshot pseudotime on 7,952 neoplastic cells (LLS origin) placed WEBS along a transcriptional continuum (GAM p < 10⁻¹⁶; Spearman ρ = 0.301, p = 2.7 × 10⁻¹⁶⁶; Fig. 5). Because the LLS origin is an analytical convention and the trend attenuated at the patient level (per-patient median ρ = 0.09, p = 0.22), we interpret this as supportive evidence for a continuous axis of state-associated variation rather than proof of a directional LLS→LDTS→LHS transition.",
"Slingshot pseudotime on 7,952 neoplastic cells placed WEBS along a transcriptional continuum (GAM p < 10⁻¹⁶; ρ = 0.301; Fig. 5). Because the trend attenuated at the patient level (p = 0.22), we interpret it as supportive evidence for continuous state-associated variation rather than proof of directional transition.")

# 3.4 (compress the rigor enumeration; keep key numbers)
s = s.replace(
"Protein-level WEBS stratified 99 CPTAC-GBM tumors. Differential abundance (limma, LHS vs LLS) revealed 107 of 111 quantified spliceosome proteins significantly up-regulated in LHS (FDR < 0.05) with none significantly down-regulated (enrichment p = 3.4 × 10⁻³⁷); notably, 110 of the 111 proteins showed positive log fold-changes, indicating a broadly shifted abundance pattern across the complex. The direction-uniform up-regulation was robust to stratification threshold (109–111/111 across 20%/25%/30% quantile cutoffs and median split; Supplementary Table S_threshold_sensitivity). Spliceosome abundance tracked protein-level WEBS (ρ = 0.79, 95% bootstrap CI 0.68–0.87) but not mRNA-level WEBS (ρ = −0.03, 95% CI −0.22 to 0.18); the layer discordance (Δρ = 0.81) was formally significant (bootstrap 95% CI 0.60–1.03, p < 0.0001), was not driven by any single protein (leave-one-protein-out ρ range 0.788–0.795), held at the per-gene level (107/109 genes showed stronger protein- than mRNA-layer correlation; Fig. S22), and persisted after controlling for proliferation and immune-content protein scores (partial ρ = 0.65, p = 9.7 × 10⁻¹³) as well as mesenchymal, hypoxia, and OPC-like programs (joint partial ρ = 0.66, p = 5.8 × 10⁻¹³). Notably, an expression-matched permutation control (10,000 random gene sets matched for proteome abundance) showed that partial protein–mRNA discordance is a generic property of abundant protein complexes (spliceosome empirical p = 0.22 by this stricter criterion; Fig. S7); the spliceosome-specific claims in this study therefore rest on the uniform direction of regulation (107/111), its enrichment statistics, and the functional layers below, rather than on discordance magnitude alone. Protein WEBS also correlated with overall histone acetylation (r = 0.62, p = 5.3 × 10⁻¹²).",
"Protein-level WEBS stratified 99 CPTAC-GBM tumors. Differential abundance revealed 107 of 111 quantified spliceosome proteins significantly up-regulated in LHS (FDR < 0.05; none down; enrichment p = 3.4 × 10⁻³⁷), with 110/111 positive — robust across stratification thresholds (109–111/111; Table S_threshold_sensitivity). Spliceosome abundance tracked protein-level WEBS (ρ = 0.79, bootstrap 95% CI 0.68–0.87) but not mRNA-level WEBS (ρ = −0.03; Δρ = 0.81, 95% CI 0.60–1.03), robust to leave-one-protein-out and to proliferation/immune/MES/hypoxia/OPC adjustment (joint partial ρ = 0.66). Notably, expression-matched permutation (10,000 random gene sets) showed partial discordance to be generic to abundant complexes (empirical p = 0.22; Fig. S7); spliceosome-specific claims therefore rest on the direction-uniform up-regulation, enrichment statistics, and the functional layers below. Protein WEBS also correlated with overall histone acetylation (r = 0.62).")

# 3.9 signature section (tighten first paragraph)
s = s.replace(
"The six-gene signature (HIF1A, SIRT2, MMP2, ITGA5, HDAC1, SPP1; fixed coefficients) stratified overall survival in TCGA-GBM (log-rank p = 0.0028), CGGA-693 (p = 3.4 × 10⁻¹⁹), and CGGA-325 (p = 7.7 × 10⁻¹⁸), and was further validated in the GLASS cohort (107 primary tumors: Cox HR = 2.14, p = 2.0 × 10⁻⁶; log-rank p = 4.2 × 10⁻⁵; Fig. 6). Two additional cohorts showed concordant directions without significance: GSE7696 (n = 84; HR = 2.59, 95% CI 0.92–7.31, p = 0.073) and GSE83300 (n = 50; HR = 2.76, 95% CI 0.01–1043, p = 0.74), the latter illustrating the uncertainty inherent to small cohorts (Fig. S11). Attenuation in IDH-wild-type-only subsets (CGGA-693 p = 0.13; CGGA-325 p = 0.62) is reported explicitly as a boundary (Fig. 6a–d).",
"The six-gene signature (HIF1A, SIRT2, MMP2, ITGA5, HDAC1, SPP1; fixed coefficients) stratified survival in TCGA-GBM (log-rank p = 0.0028), CGGA-693 (p = 3.4 × 10⁻¹⁹), CGGA-325 (p = 7.7 × 10⁻¹⁸), and GLASS (HR = 2.14, p = 2.0 × 10⁻⁶; Fig. 6a–d). Two additional cohorts showed concordant directions without significance (GSE7696 p = 0.073; GSE83300 p = 0.74; Fig. S11), and attenuation in IDH-wild-type subsets is reported explicitly as a boundary.")

# ---------- 6. compression: Discussion 4.1 ----------
s = s.replace(
"This study identifies a writer-dominant expression state associated with post-transcriptional spliceosome accumulation and altered splicing output in GBM, accompanied by coordinated changes in post-translational modification. The evidence chain is layered: (i) spliceosome proteins accumulate without proportional mRNA induction (ρ = 0.79 vs −0.03); (ii) their lysine acetylation is enriched (OR = 7.36 raw; 2.29 after abundance correction), shows a significant ordered trend across strata (43/43 sites), and couples to abundance (median ρ = 0.62); (iii) alternative splicing output is systematically altered (4,387 events), enriched for splicing-related genes; (iv) spliceosome phosphorylation increases with abundance while occupancy-level phosphorylation decreases, alongside elevated splicing-kinase abundance; and (v) connectivity-based sensitivity analysis places spliceosome genes disproportionately downstream of acetylation regulators (regulator-level p = 6.1 × 10⁻⁵). These findings are supported by complementary molecular assay layers and additional cohort-level analyses, although several omics measurements were obtained from the same underlying patient cohorts.",
"This study identifies a writer-dominant expression state associated with post-transcriptional spliceosome accumulation and altered splicing output in GBM. The evidence chain is layered: protein-level accumulation without proportional mRNA induction (i); enriched, abundance-corrected acetylation coupled to protein abundance (ii); systematically altered splicing output enriched for splicing-related genes (iii); abundance-dependent phosphorylation changes with elevated splicing-kinase abundance (iv); and connectivity-based sensitivity of spliceosome genes to acetylation-regulator perturbation (v). Several omics layers derive from the same underlying cohorts, which we account for in our interpretation.")

# ---------- 7. Additional Information: BJC order + Consent + ORCID ----------
addinfo = """# Acknowledgements

The authors thank the data-generating consortia (TCGA, CPTAC, CGGA, GLASS, DepMap, GDSC) and the original study authors.

# Authors' contributions

**Jia Chen**: Conceptualization, Data curation, Formal analysis, Investigation, Methodology, Software, Validation, Visualization, Writing – original draft. **Junhan Yang**: Validation, Writing – review & editing. **Bing Yang**: Writing – review & editing. **Tianyi Zhou**: Writing – review & editing. **Peijun Zhao**: Writing – review & editing.

# Ethics approval and consent to participate

This study involved only secondary analyses of publicly available, de-identified datasets (TCGA, CPTAC, CGGA, GEO, GLASS, DepMap, GDSC2, TCGASpliceSeq). No new human participants, tissue samples, or identifiable data were collected. All datasets were accessed and analyzed in accordance with their respective data-use agreements.

# Consent for publication

Not applicable.

# Data availability

All datasets analyzed in this study are available from their respective repositories, subject to the access requirements and terms of use established by the original data providers: GSE138794, GSE194329, GSE131928, GSE7696, GSE83300 (GEO, open access); TCGA cohorts (GDC and the UCSC Xena GDC hub, STAR-TPM release); CGGA-693/325 (cgga.org.cn, registration and data-use agreement required); CPTAC-GBM (LinkedOmics); TCGASpliceSeq (bioinformatics.mdanderson.org); DepMap 24Q2 (depmap.org); GDSC2 (cancerrxgene.org). GLASS data were accessed via Synapse (syn17038081) under the consortium's data-use terms and are not redistributed. Analysis code and the WEBS implementation are publicly available at https://github.com/huozheshuashuai/GBM_WEBS; the submission version corresponds to the tagged release indicated in the repository README.

# Competing interests

The authors declare that they have no competing interests.

# Funding information

The authors received no specific funding for this work.

# References"""

s = re.sub(r"# Ethics approval and consent to participate.*?# References", addinfo, s, flags = re.S)

# ---------- 8. AI disclosure in Methods ----------
s = s.replace("R 4.5.2; full pipeline (50+ scripts) and the WEBS standalone implementation are available in the project repository.",
              "R 4.5.2; full pipeline (60+ scripts) and the WEBS standalone implementation are available in the project repository. During preparation, generative AI tools (large language models) were used for language editing, code formatting, and pipeline reconstruction assistance; all analyses were executed by the authors, and all code, outputs, and scientific conclusions were independently verified by the authors against the underlying data.")

# ---------- 9. Figure 6/7 merge: replace legends ----------
new_legends = """# Figure Legends

**Figure 1. Study design and analytical workflow.** Data sources, WEBS computation, and the analysis tiers (discovery, functional characterization, clinical association, candidate vulnerability nomination).

**Figure 2. The Writer–Eraser Balance Score (WEBS) landscape.** (a) Cell-type annotation of 35,823 GBM cells. (b) WEBS projected onto the UMAP. (c) WEBS by cell type. (d) WEBS by IQR group (LLS/LDTS/LHS; Kruskal–Wallis χ² = 30,439, df = 2, p < 2.2 × 10⁻¹⁶ [reportable lower bound]; point shapes distinguish groups for grayscale readability).

**Figure 3. Protein-level spliceosome accumulation is associated with the writer-dominant state (CPTAC-GBM).** (a) Proteome differential abundance (high- vs low-WEBS); spliceosome up-regulated proteins highlighted (triangles; 107/111, FDR < 0.05). (b) Protein-level WEBS vs spliceosome protein score (Spearman ρ = 0.79, 95% CI 0.68–0.87). (c) Per-gene mRNA- vs protein-layer correlation with WEBS (107/109 above diagonal). (d) Robustness across stratification thresholds (109–111/111 up-regulated).

**Figure 4. Coordinated post-translational modification of the spliceosome.** (a) Top differential spliceosome acetylation sites, samples ordered by WEBS (enrichment OR = 7.36; abundance-corrected OR = 2.29). (b) Phosphosite volcano (high- vs low-WEBS): 198/319 spliceosome sites significant (197 increased); after parent-protein correction, significant sites predominantly decreased in abundance-adjusted signal (118/124). (c) Splicing-kinase abundance vs WEBS. (d) GENIE3-based connectivity sensitivity analysis: relative absolute predicted shift of spliceosome vs other genes per regulator perturbation (dashed line = ratio 1; purple = writers, teal = erasers).

**Figure 5. Altered alternative splicing output.** (a) 28,071 AS events vs bulk WEBS (TCGA-GBM); 4,387 significant (dashed lines = |ρ| = 0.2 threshold). (b) Event-type composition. (c) GO enrichment of genes harboring WEBS-associated splicing events. (d) Cross-cohort replication in TCGA-LGG (83% sign-concordant of 4,285 testable events; effect-size ρ = 0.70).

**Figure 6. Clinical validation and longitudinal association.** (a–d) Kaplan–Meier curves with number-at-risk tables for the six-gene signature (fixed coefficients; cohort-median split; log-rank p in panel, HR with 95% CI from group Cox models): TCGA-GBM (training), CGGA-693, CGGA-325, GLASS primary tumors. (e) Primary-tumor WEBS vs observed interval to first recurrence (GLASS, first record per patient, n = 173; exploratory). (f) WEBS in paired primary–recurrence samples (179 pairs from 148 patients).

# Tables"""
s = re.sub(r"# Figure Legends\n\n.*?# Tables", new_legends, s, flags = re.S)

# ---------- 10. figure reference remap: Fig. 6d -> Fig. 6; drop Fig. 7 refs ----------
s = s.replace("Fig. 7a", "Fig. 6e").replace("Fig. 7b", "Fig. 6f")
s = s.replace("(Fig. 7c,d)", "(Fig. S35)")
s = s.replace("Fig. 14a", "Fig. 6d")
# Fig 7 panel references inside text
s = s.replace("Top WEBS-sensitizing single agents", "Top WEBS-sensitizing single agents")

# ---------- 11. remove the Figures section (BJC: figures uploaded separately) ----------
s = re.sub(r"# Figures\n\n.*", "", s, flags = re.S)

# ---------- 12. corresponding author ORCID placeholder ----------
s = s.replace("\\* Corresponding author: Jia Chen, Department of Bioinformatics, Chengdu Neusoft University, Chengdu, China. Email: 1176625889@qq.com",
              "\\* Corresponding author: Jia Chen, Department of Bioinformatics, Chengdu Neusoft University, Chengdu, China. Email: 1176625889@qq.com. ORCID: [TO BE ADDED BY AUTHOR]")

open("MANUSCRIPT_BJC.md", "w", encoding="utf-8").write(s)

# word count of body sections
def wc(t):
    return len(re.findall(r"[\w’'-]+", t))
parts = re.split(r"(?m)^(##? .+)$", s)
total = 0
for i, p in enumerate(parts):
    if p.startswith("# 1.") or p.startswith("# 2.") or p.startswith("# 3.") or p.startswith("# 4.") or p.startswith("# 5.") or p.startswith("## "):
        pass
body = re.search(r"# 1. Introduction(.*?)# References", s, re.S).group(1)
print("body words:", wc(body))
print("MANUSCRIPT_BJC.md written")
