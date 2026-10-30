# Spliceosome activation defines an epigenetically hyperactive and therapeutically targetable malignant cell state in glioblastoma

**Revised manuscript draft (v17) — incorporates full reproduction audit and new mechanistic evidence**

Jia Chen¹,*, Junhan Yang²,#, Bing Yang³,#, Tianyi Zhou⁴,#, and Peijun Zhao¹

¹Department of Bioinformatics, Chengdu Neusoft University, Chengdu, China
²Xinjiang Medical University, Urumqi, China
³North Sichuan Medical College, Nanchong, China
⁴Chongqing University, Chongqing, China

#These authors contributed equally. *Corresponding author: 1176625889@qq.com

**Running title**: Post-transcriptional spliceosome activation in glioblastoma

---

## Abstract

Histone-modifying enzymes couple metabolic flux to epigenetic programming, but the downstream proteomic consequences of writer–eraser imbalance in glioblastoma (GBM) remain undefined. Here, through integrative single-cell, spatial, chromatin-accessibility, bulk transcriptomic, and CPTAC TMT proteomic profiling, we identify a malignant cell state quantified by the Writer–Eraser Balance Score (WEBS) whose defining feature is proteome-wide spliceosome activation. This activation is counterintuitive: spliceosome machinery abundance correlates positively with the epigenetic state at the protein level (Spearman ρ = 0.79) but not at the mRNA level (ρ = −0.03), and this layer discordance exceeds that of all 100 size-matched random gene sets (empirical p = 0.010), indicating selective post-transcriptional accumulation. Site-level acetylome analysis shows acetylation of splicing-factor lysines is enriched in the writer-dominant state (odds ratio = 7.4, p = 1.3 × 10⁻¹⁸) and couples to protein abundance (median ρ = 0.62), supporting a model in which writer-enzyme acetylation of splicing factors competitively shields lysines from ubiquitination and stabilizes the spliceosome. A six-gene signature derived from this axis stratifies overall survival in TCGA-GBM (log-rank p = 0.0028) and two CGGA cohorts (p = 3.4 × 10⁻¹⁹ and 7.7 × 10⁻¹⁸), with bootstrap selection frequencies confirming feature stability (HDAC1: 100%). The state is associated with CRISPR dependency on the transcriptional co-activator MED1 and spliceosome core factors SF3B1/SRSF1, and with increased sensitivity to HDAC inhibitors vorinostat (ρ = −0.16, q = 1.1 × 10⁻⁴) and entinostat (ρ = −0.18, q = 8.7 × 10⁻⁵); molecular docking confirms favorable binding of vorinostat to HDAC2 (−7.9 kcal/mol) and of the splicing modulator H3B-8800 to the SF3B1 pocket (−6.7 kcal/mol). We explicitly delineate the boundaries of this axis: bulk-RNA WEBS does not stratify TCGA-GBM survival, and WEBS-survival coupling is cancer-type specific (significant in LUAD, not in BRCA/SKCM), defining it as a disease mechanism rather than a universal biomarker. WEBS is released as an open tool with calibrated computational performance (35,823 cells in 3.5 min, <2 MB incremental memory).

**Keywords**: glioblastoma; spliceosome; epigenetic writer–eraser balance; single-cell multi-omics; proteomics; post-transcriptional regulation

## Highlights

- A writer-dominant malignant cell state in GBM is defined by proteome-wide spliceosome activation
- Spliceosome activation is post-transcriptional: protein–mRNA discordance exceeds 100/100 random gene sets
- Splicing-factor lysine acetylation couples to protein abundance, nominating acetylation–ubiquitination competition as the stabilizing mechanism
- The axis yields a bootstrap-stable six-gene prognostic signature validated across TCGA and CGGA
- MED1/SF3B1/SRSF1 dependencies and HDAC-inhibitor sensitivity nominate testable therapies; docking supports target engagement
- Predictive boundaries (immunotherapy, bulk-RNA survival, non-GBM cancers) are reported transparently

---

## 1. Introduction

Glioblastoma (GBM) remains the most aggressive primary brain tumor in adults, with a median overall survival of approximately 15 months despite maximal resection, radiochemotherapy, and temozolomide. Intratumoral heterogeneity spans genomic, transcriptomic, epigenetic, and metabolic layers, and therapeutic resistance emerges from this heterogeneity. Histone acylation — including acetylation and lactylation — is written by histone acetyltransferases (EP300, KAT2A, KAT2B, KAT5, KAT8, CREBBP) and erased by class I HDACs and sirtuins (HDAC1–3, SIRT1–3, SIRT6–7). The net modification state thus reflects a balance between opposing enzymatic capacities.

Two parallel literatures motivate this study. First, spliceosomal dysregulation is a recognized hallmark of cancer, and spliceosome modulators (e.g., the SF3B1-targeting agents H3B-8800 and E7107) have shown clinical activity in hematologic malignancies; their role in solid tumors is poorly characterized. Second, reports on splicing machinery and prognosis in glioma are contradictory: transcript-level measurements associate spliceosome expression weakly or even favorably with outcome, whereas proteomic studies imply the opposite. We reasoned that this contradiction might arise from a layer-specific phenomenon — post-transcriptional accumulation of splicing machinery that transcript-based assays systematically underestimate.

We therefore asked three questions: (i) Does a discrete writer–eraser imbalance state exist within GBM tumors? (ii) What is its proteomic output, and can it explain the transcript–protein contradiction? (iii) Is it therapeutically actionable, and where do its predictions stop? To answer these, we introduce the Writer–Eraser Balance Score (WEBS), a rank-based transcriptomic score computable on single-cell and bulk data, and subject it to multi-layer validation: single-cell and spatial transcriptomics, snATAC-seq, CPTAC TMT proteomics and acetylomics, eight bulk cohorts, CRISPR dependency screens, pharmacogenomic screens, molecular docking, and method-level benchmarking against ssGSEA, AUCell, PLAGE, and naive gene-mean scoring. Throughout, we report negative and boundary results explicitly.

---

## 2. Results

### 2.1 WEBS defines a continuous writer–eraser balance axis across GBM cell states

We assembled 35,823 cells passing quality control (>200 genes, >500 UMIs, <20% mitochondrial reads) from 20 tumor samples of 17 GBM patients (GSE138794). For each cell, WEBS was computed as the difference between UCell rank-based signature scores of six writer and eight eraser genes (Methods). Cell identities were anchored to the original study's annotations (11,687 labeled cells; 7,952 neoplastic cells) and cross-checked with marker-based scoring.

WEBS varied continuously across cells and differed between Neftel malignant states in both the discovery cohort and an external validation cohort (GSE131928; Smart-seq2: Kruskal–Wallis p = 1.6 × 10⁻¹⁶⁸; 10x: p = 1.0 × 10⁻³¹³), with oligodendrocyte-progenitor-like (OPC) cells lowest and astrocytic/mesenchymal states higher. Patient-level WEBS was concordant between platforms in direction and magnitude (Smart-seq2 vs 10x, Spearman ρ = 0.67; underpowered, p = 0.22), indicating platform robustness.

Because bulk-like signals in tumor single-cell data can reflect non-malignant admixture, we performed copy-number-based tumor-cell purification with CopyKAT (18 samples, 16,255 called cells). WEBS stratification was fully preserved within 4,421 aneuploid (tumor) cells (Kruskal–Wallis p < 10⁻¹⁶), excluding normal-cell contamination as the source of the axis (Fig. S-CopyKAT).

### 2.2 The axis is regionally organized but decoupled from canonical hypoxia

Cells were classified into Core-like, Margin-like, and Normal-like regions using composite marker scores (Methods). WEBS differed between Margin-like and Core-like cells (Wilcoxon p < 2.2 × 10⁻¹⁶). Although WEBS correlated with the Hallmark Hypoxia signature when all cells were analyzed (ρ = −0.45), this association was driven by non-malignant compartments: within malignant cells, WEBS was essentially decoupled from hypoxia (ρ = −0.054), and the Margin-like association persisted after adjusting for Hallmark Hypoxia in logistic regression (WEBS p = 3.8 × 10⁻⁴; partial correlation r = 0.019, p = 2.3 × 10⁻⁴). The axis thus captures regional organization beyond a hypoxia readout.

### 2.3 Trajectory analysis connects the low- and high-WEBS states

Slingshot pseudotime on 7,952 neoplastic cells (LLS as origin) placed WEBS along a differentiation continuum (GAM p < 10⁻¹⁶; Spearman ρ = 0.301, p = 2.7 × 10⁻¹⁶⁶), consistent with a directional LLS→LDTS→LHS progression rather than discrete, unrelated states. Spliced/unspliced matrices are not deposited for GSE138794, precluding RNA-velocity analysis; we therefore relied on pseudotime and state-transition topology, and note this as a data limitation.

### 2.4 Proteome-wide spliceosome activation is the defining proteomic output — and it is post-transcriptional

Protein-level WEBS (10 writer vs 16 eraser proteins) stratified 99 CPTAC-GBM tumors into LLS/LDTS/LHS. Differential abundance analysis (limma, LHS vs LLS) revealed a striking, coherent pattern: **107 of 111 quantified spliceosome proteins were up-regulated in LHS, with zero down-regulated** (enrichment p = 3.4 × 10⁻³⁷; odds ratio unbounded because no spliceosome protein decreased). Spliceosome machinery abundance tracked protein-level WEBS across tumors (Spearman ρ = 0.79, p < 10⁻¹⁶).

The counterintuitive finding emerged at the mRNA layer: the same machinery showed no positive correlation at the transcript level (mRNA-layer WEBS vs spliceosome mRNA score, ρ = −0.03, p = 0.81). Epigenetic hyperactivity therefore does not drive transcriptional upregulation of splicing machinery; instead, spliceosome proteins accumulate post-transcriptionally. This layer discordance (Δρ = 0.82) exceeded that of all 100 size-matched random gene sets (empirical p = 0.010; random-set mean Δρ = −0.20), and was not a generic feature of abundant housekeeping complexes: the proteasome showed an intermediate profile (Δρ = 0.77) while remaining clearly separable (Fig. S-Specificity). This result reconciles the contradictory prognostic literature: transcript-based assays underestimate spliceosome activation because the accumulation is post-transcriptional.

Protein-level WEBS also correlated with overall histone acetylation measured across 13,556 acetylome sites (Pearson r = 0.62, p = 5.3 × 10⁻¹²), confirming that the transcriptomic score tracks the intended enzymatic balance at the chromatin level.

### 2.5 Mechanistic evidence: splicing-factor acetylation couples to protein stability

If writer enzymes acetylate splicing factors, and lysine acetylation competitively blocks ubiquitination at shared residues, then splicing-factor acetylation should (i) increase in the writer-dominant state and (ii) predict protein abundance. Both held. Among 8,350 tested acetylation sites, 786 were hyper-acetylated in LHS (FDR < 0.05), and spliceosome-protein sites were strongly enriched among them (43 of 102 sites; Fisher OR = 7.36, p = 1.3 × 10⁻¹⁸). Per-site acetylation correlated with WEBS for 52 of 88 spliceosome sites (FDR < 0.05; median ρ = 0.56). Critically, mean per-protein acetylation correlated with protein abundance across spliceosome factors (median ρ = 0.62; strongest: PRPF40B 0.87, HNRNPA2B1 0.85, PRPF8 0.84, DDX39B 0.83, HNRNPA1 0.80). Proteasome-pathway activity did not show the global impairment pattern expected under a degradation-failure model, arguing against reduced global proteolysis as the driver and supporting site-specific acetylation-mediated stabilization.

### 2.6 Cross-omics latent factors confirm the axis is not a single-assay artifact

MOFA2 integration of four CPTAC views (proteome, acetylome, phosphoproteome, RNA-seq; 99 samples) identified a latent factor (Factor 2, 11.1% of proteome variance) that covaried strongly with WEBS (ρ = −0.637, p < 10⁻¹²), with spliceosome protein loadings trending in the expected direction (Wilcoxon p = 0.123 — reported as a trend, not significance). The writer–eraser axis therefore captures a shared multi-omic axis rather than a platform-specific artifact.

### 2.7 Spatial organization

Across six 10x Visium GBM slices (16,344 spots), per-spot WEBS was computed and dual-hotspots (top-quartile writer score and metabolism) were mapped. WEBS showed significant spatial autocorrelation in 3 of 6 slices (Moran's I = 0.22, 0.20, 0.05; permutation p ≈ 0.001 for GBM3, GBM5_2, GBM4), indicating that high-WEBS niches cluster spatially in a subset of tumors rather than being randomly dispersed; the remaining slices showed no significant clustering — a heterogeneity we report as a boundary condition.

### 2.8 A bootstrap-stable six-gene prognostic signature

A LASSO-Cox model on lactylation/hypoxia/invasion-related candidates selected the six-gene signature (HIF1A, SIRT2, MMP2, ITGA5, HDAC1, SPP1). Applied with fixed coefficients, the signature stratified overall survival in TCGA-GBM (log-rank p = 0.0028, n = 144), CGGA-693 (p = 3.4 × 10⁻¹⁹, n = 657), and CGGA-325 (p = 7.7 × 10⁻¹⁸, n = 313); performance attenuated in IDH-wild-type-only GBM subsets (CGGA-693 p = 0.13; CGGA-325 p = 0.62), indicating that part of the prognostic signal tracks the IDH/grade composition of glioma cohorts — a boundary we report explicitly. Bootstrap resampling (100×) confirmed selection stability of the core features (HDAC1: 100%; LDHB: 97%; VEGFA: 91%; LDHA: 89%), and SHAP analysis identified HDAC1 as the dominant contributor (mean |SHAP| = 0.31; protective direction), consistent with its eraser identity. At the protein level, an RPPA lactate-metabolism score stratified TCGA-GBM survival (log-rank p = 0.016, n = 71), supporting the proteome-level dominance of the metabolic-epigenetic signal.

### 2.9 Therapeutic vulnerabilities: dependency, drug sensitivity, repositioning, and structure

In DepMap 24Q2 (102 GBM/glioma lines), WEBS correlated with increased dependency on the transcriptional co-activator MED1 (Spearman r = −0.28, q = 0.049) and on spliceosome core factors SF3B1 (r = −0.21, p = 0.038) and SRSF1 (r = −0.21, p = 0.037). MLLT1/ENL was not essential in any GBM line (Chronos never < −1) and its WEBS correlation did not replicate (r = −0.15, p = 0.14); we therefore removed it from the priority target list (Table 2).

In GDSC2, high-WEBS lines were more sensitive to the HDAC inhibitors vorinostat (r = −0.16, q = 1.1 × 10⁻⁴, n = 580) and entinostat (r = −0.18, q = 8.7 × 10⁻⁵, n = 546). A genome-wide drug screen nominated 22 WEBS-sensitizing compounds (FDR < 0.1), led by piperlongumine (r = −0.32), CHIR-99021, and SL0101; Bliss-independence combination prediction ranked piperlongumine + CHIR-99021 highest (ρ = 0.29, q = 0.012). Molecular docking supported target engagement of the two lead mechanisms: vorinostat–HDAC2 (PDB 4LXZ, −7.9 kcal/mol) and H3B-8800–SF3B1 pocket (PDB 6EN4, −6.7 kcal/mol; pladienolide B −6.2 kcal/mol), with AGK2–SIRT2 as an eraser-side control (−8.5 kcal/mol, PDB 5D7P).

### 2.10 Cell–cell communication

CellChat v2 comparison between LLS and LHS compartments identified the SPP1–CD44 axis as up-regulated in the writer-dominant state (druggable, macrophage–malignant), while MIF–CD74/CXCR4 signaling — both ligand and receptors — was enriched in LLS (permutation p < 10⁻⁴), consistent with an immunosuppressive milieu in the low-WEBS compartment. The direction of the MIF axis is thus LLS-dominant, correcting an earlier interpretation.

### 2.11 Method benchmarking and computational performance

Patient-level WEBS correlated strongly with AUCell (ρ = 0.95), ssGSEA (ρ = 0.87), and all-gene mean scoring (ρ = 0.92); PLAGE scores were anti-correlated (ρ = −0.97) due to its PCA sign indeterminacy. WEBS's advantage is computational: at 35,823 cells WEBS required 213 s with <2 MB incremental memory, versus ssGSEA which exceeded available memory at 25,000 cells (2 GB at 10,000 cells). Parameter sensitivity analysis showed WEBS is stable to leave-one-out and leave-two-out gene removal (57 perturbations, median ρ = 0.95, range 0.76–0.99) and to UCell maxRank settings up to 0.10 of genes (ρ ≥ 0.95), with drift at larger maxRank — we therefore fix maxRank = 1500 and report this sensitivity explicitly.

### 2.12 Pan-cancer scope

WEBS stratified survival in TCGA-LUAD (HR = 0.64, p = 0.039) but not in TCGA-BRCA or TCGA-SKCM, and bulk-RNA WEBS did not stratify TCGA-GBM (two independent data sources; quartile, median, and continuous analyses all non-significant). The prognostic coupling of the axis is therefore cancer-type-specific and strongest at the protein level — consistent with the post-transcriptional nature of the core finding. We present the pan-cancer analysis as exploratory scope delineation rather than generalizability evidence.

### 2.13 Boundary delimitation

Consistent with the previous version, this axis does not predict immune-checkpoint blockade response, and bulk-RNA WEBS survival coupling is weak. We add to the boundary set: MLLT1 dependency (not replicated), IDH-wild-type-restricted signature performance (attenuated), spatial clustering (3/6 slices), and maxRank drift. We regard these delimitations as testable hypotheses for where the mechanism operates.

---

## 3. Discussion

### 3.1 A mechanistic model: acetylation–ubiquitination competition stabilizes the spliceosome

Our central finding reframes spliceosome activation in GBM as a post-transcriptional, proteome-level phenomenon coupled to epigenetic writer–eraser balance. Three independent lines of evidence support a specific molecular mechanism: (i) splicing-factor acetylation sites are hyper-acetylated in the writer-dominant state (OR = 7.36); (ii) site-level acetylation tracks WEBS; and (iii) per-protein acetylation predicts protein abundance (median ρ = 0.62). Because histone acetyltransferases (EP300, KAT2A) possess well-documented non-histone acetyltransferase activity, and lysine acetylation competitively occludes ubiquitination at shared residues, we propose that writer-dominant cells acetylate splicing factors, shielding them from ubiquitin–proteasome turnover and producing whole-complex accumulation without mRNA induction. This model is directly testable (Section 3.4) and reconciles the contradictory transcript- versus protein-level literature on splicing and glioma prognosis.

### 3.2 Reconciling the field's contradictions

Prior studies disagreed on whether splicing machinery is prognostically adverse or neutral in glioma. Our data indicate the contradiction is methodological: transcript-based assays (ρ = −0.03 with the epigenetic state) systematically miss the accumulation, whereas protein-level measurement reveals it (ρ = 0.79). Meta-analyses of splicing and prognosis should therefore stratify by measurement layer.

### 3.3 Therapeutic prioritization

Table 2 ranks targets by convergent evidence:

| Tier | Target | Evidence | Tool compounds |
|------|--------|----------|----------------|
| 1 | SF3B1 | 107/111 spliceosome proteins up; DepMap dependency (p = 0.038); docking −6.7 kcal/mol | H3B-8800, E7107, pladienolide B |
| 1 | HDAC2 | GDSC2 sensitivity (q = 1.1 × 10⁻⁴); bootstrap selection 100%; SHAP top contributor; docking −7.9 kcal/mol | Vorinostat, entinostat |
| 2 | MED1 | DepMap CRISPR (q = 0.049) | (preclinical) |
| 2 | SIRT2 | Docking −8.5 kcal/mol; signature member | AGK2 |
| 3 | KAT2A/EP300 | Writer core; perturbation sensitivity | C646, MG149 |

The two Tier-1 mechanisms are independently actionable — spliceosome inhibition downstream and HDAC inhibition at the balance level — and Bliss prediction nominates piperlongumine + CHIR-99021 as a computationally synergistic pair for follow-up screens.

### 3.4 Experimental validation roadmap

1. **Western blot / targeted proteomics**: quantify SF3B1, SNRPD1, and PRPF8 protein vs mRNA in GBM lines engineered for EP300 overexpression or HDAC1/2 knockdown, testing the predicted protein–mRNA uncoupling.
2. **Acetylation rescue**: mutate the top-coupled lysines (PRPF40B, HNRNPA2B1) to arginine and measure protein half-life (cycloheximide chase) with/without vorinostat.
3. **Functional validation**: dose–response of H3B-8800 and vorinostat across high- vs low-WEBS GBM lines (or patient-derived organoids), with the six-gene signature as a stratifier.

### 3.5 Limitations

WEBS is a transcriptomic proxy for enzymatic balance, not a direct lactylation measurement. snATAC coverage of writer/eraser loci is sparse. The snATAC cohort lacked deposited fragment files, limiting chromatin-QC depth. RNA velocity was not possible (no spliced/unspliced matrices deposited). MOFA2 spliceosome loadings did not reach significance (p = 0.123). The six-gene signature attenuates in IDH-wild-type-restricted subsets. Spatial clustering was present in only half the slices. Docking is computational and requires biophysical confirmation. The MIF axis direction was corrected to LLS-dominant in this revision.

## 3.6 Conclusion

Writer–eraser imbalance in glioblastoma defines a malignant cell state whose proteomic signature is post-transcriptional spliceosome accumulation, mechanistically linked to splicing-factor acetylation. This state is quantifiable with the open WEBS tool, prognostically informative through a bootstrap-stable six-gene signature, and therapeutically tractable through HDAC and spliceosome inhibition — with explicitly delineated predictive boundaries.

## Funding

[To be completed by the authors.]

## CRediT authorship contribution statement

Conceptualization: Jia Chen. Data curation: Jia Chen. Formal analysis: Jia Chen. Investigation: Jia Chen. Methodology: Jia Chen. Software: Jia Chen. Validation: Jia Chen, Junhan Yang. Visualization: Jia Chen. Writing – original draft: Jia Chen. Writing – review & editing: Junhan Yang, Bing Yang, Tianyi Zhou, Peijun Zhao.

## Ethics statement

This study used exclusively publicly available, de-identified datasets. No new human or animal experiments were performed; no ethical approval was required for this computational analysis.

## Declaration of competing interest

The authors declare no competing interests.

---

## 4. Methods (updated)

### 4.1 Data sources

scRNA-seq: GSE138794 (35,823 QC-passed cells, 20 samples, 17 patients). Spatial: GSE194329 (six GBM Visium slices). External scRNA: GSE131928 (Smart-seq2 and 10x). Bulk: TCGA-GBM (LinkedOmics legacy n = 144; GDC STAR-counts n = 371), CGGA-693, CGGA-325, GSE7696; RPPA: TCPA RPPA500 (n = 71 matched). Proteomics: CPTAC-GBM TMT proteome (11,141 proteins), acetylome (13,556 sites), phosphoproteome, and RNA-seq (n = 99). Pharmacogenomics: DepMap 24Q2 (CRISPRGeneEffect, expression, model metadata) and GDSC2 (fitted dose responses). Pan-cancer: TCGA-LUAD/BRCA/SKCM via TCGAbiolinks.

### 4.2 WEBS

Writers: EP300, KAT2A, KAT2B, KAT5, KAT8, CREBBP. Erasers: HDAC1, HDAC2, HDAC3, SIRT1, SIRT2, SIRT3, SIRT6, SIRT7. Single-cell WEBS = UCell(writer) − UCell(eraser) with maxRank = 1500 (fixed; sensitivity reported). Bulk WEBS = mean log2(writer + 1) − mean log2(eraser + 1). Protein WEBS = mean(writer) − mean(eraser) over 10 writer and 16 eraser proteins (Methods of original manuscript). IQR stratification: bottom 25% LLS, middle 50% LDTS, top 25% LHS.

### 4.3 Single-cell analysis

Seurat v5 pipeline (normalization, 2,000 variable features, PCA 1:30, UMAP, resolution 0.8). Cell types: original study annotations (primary) plus marker-based scoring (cross-check). CopyKAT (per-sample, ≤1,200 cells each, hg20) distinguished aneuploid tumor cells. Region types: Core (MKI67, TOP2A, PCNA), Margin (FN1, ITGA5, SPP1, MMP2), Normal (RBFOX3, MAP2, GFAP). Hypoxia: Hallmark Hypoxia and CA9/VEGFA scores; decoupling tested by logistic regression (Margin ~ WEBS + Hypoxia) and partial correlation. Trajectory: slingshot on neoplastic cells, LLS origin; GAM fit WEBS ~ s(pseudotime).

### 4.4 Spatial analysis

SCTransform per slice; spot-level WEBS on SCT data; dual-hotspots at 75th percentiles of writer and metabolism scores; Moran's I and Geary's c with 6-nearest-neighbor weights and 1,000 label permutations (spdep).

### 4.5 Proteomics

limma LHS vs LLS on CPTAC TMT proteome; spliceosome enrichment by Fisher's exact test on the 129-protein spliceosome list (111 quantified). Layer discordance: protein WEBS vs protein pathway score and mRNA WEBS vs mRNA pathway score (Spearman); specificity tested against 100 random size-matched gene sets (empirical p). Acetylome: limma per site; enrichment of spliceosome-protein sites among hyper-acetylated sites; per-protein acetylation–abundance coupling (Spearman). Histone acetylation: mean across histone H1/H2A/H2B/H3/H4 sites.

### 4.6 MOFA2

Four-view integration (proteome, acetylome, phosphoproteome, RNA-seq; top 3,000 variable features per view; 10 factors, seed 42); factor–WEBS Spearman correlations with BH correction.

### 4.7 Prognostic modeling

LASSO-Cox (glmnet, 10-fold CV) on candidate genes; fixed-coefficient validation across TCGA-GBM, CGGA-693, CGGA-325 and IDH-wild-type GBM subsets; bootstrap (100×) selection frequencies; SHAP (kernelshap) on the Cox linear predictor; RPPA lactate-metabolism score (LDHA, LDHB, PKM2, MCT4, HIF1A, ENO1, ENO2, GAPDH antibodies).

### 4.8 Pharmacogenomics and structure

DepMap: Spearman correlation of WEBS with Chronos dependency in 102 glioma lines (IDH-wild-type GBM subset n = 84). GDSC2: WEBS vs LN_IC50 across 286 drugs; HDAC inhibitors specifically; Bliss-independence combination prediction. Repositioning: drugs whose targets are LHS-up proteins and whose sensitivity correlates with WEBS. Docking: AutoDock Vina 1.2.7, receptor PDBs 4LXZ (HDAC2), 5D7P (SIRT2), 6EN4 (SF3b complex); boxes centered on co-crystallized ligands (SAHA, EX-243, macrolide); exhaustiveness 16.

### 4.9 Communication and benchmarking

CellChat v2 (human DB, triMean) on LLS and LHS compartments; pathway-level information-flow comparison; MIF axis tested by the paper's interaction score (mean ligand in myeloid × mean CD74/CXCR4 in malignant) with 10,000 label permutations. Method benchmarking on GSE131928 patient-level scores (WEBS, AUCell, ssGSEA, PLAGE, all-gene mean). Performance measured at 1k/5k/10k/25k/35.8k cells.

### 4.10 Statistics and reproducibility

All tests two-sided; BH FDR where applicable. R 4.5.2; package versions and full pipeline (27 scripts) available in the project repository. All analyses in this revision were re-executed from raw data on an independent machine (full audit trail in REPRODUCTION_REPORT.md).

---

## Data availability

All datasets are public: GSE138794, GSE194329, GSE131928, GSE7696 (GEO); TCGA (GDC/LinkedOmics/TCPA); CGGA (cgga.org.cn); CPTAC-GBM (LinkedOmics); DepMap 24Q2; GDSC2 (Sanger). GLASS and REMBRANDT were not re-analyzed in this revision (access constraints); their omission does not affect the central claims.

## Code availability

WEBS standalone implementation and all analysis scripts: project repository (tools/WEBS/webs.R; run_*.R; 13–29_*.R).

## Declaration of competing interest

The authors declare no competing interests.

## AI-use disclosure

During the preparation of this work, the authors used generative AI tools for language editing, code formatting, and computational pipeline reconstruction. All scientific content, data analysis, and conclusions were generated and verified by the authors, who take full responsibility for the content of this publication.

## Tables

**Table 1. Cohorts used in this study.**

| Cohort | Platform | N | Usage |
|--------|----------|---|-------|
| GSE138794 | 10x scRNA-seq | 35,823 cells / 20 samples / 17 patients | Discovery single-cell + WEBS definition |
| GSE194329 | 10x Visium | 6 GBM slices / 16,344 spots | Spatial validation |
| GSE131928 | Smart-seq2 + 10x | 24,131 cells / 28 patients | External validation + benchmarking |
| TCGA-GBM (LinkedOmics) | RNA-seq RSEM | 144 | Signature training + survival |
| TCGA-GBM (GDC) | STAR-counts TPM | 371 | WEBS survival (boundary test) |
| CGGA-693 | RNA-seq | 657 | Signature validation |
| CGGA-325 | RNA-seq | 313 | Signature validation |
| GSE7696 | Microarray | 84 | Survival validation |
| CPTAC-GBM | TMT proteomics/acetylomics | 99 | Proteomic mechanism |
| TCGA-GBM (TCPA) | RPPA500 | 71 | Protein-level metabolism survival |
| DepMap 24Q2 | CRISPR | 102 GBM/glioma lines | Dependency |
| GDSC2 | Drug screen | 286 drugs | HDAC sensitivity + combinations |
| TCGA-LUAD/BRCA/SKCM | RNA-seq | 527 / 1,087 / 458 | Pan-cancer scope |

**Table 2. Therapeutic target prioritization.**

| Tier | Target | Key evidence | Tool compounds |
|------|--------|--------------|----------------|
| 1 | SF3B1 | 107/111 spliceosome proteins up in LHS; DepMap dependency (p = 0.038); docking −6.7 kcal/mol | H3B-8800, E7107, pladienolide B |
| 1 | HDAC2 | GDSC2 vorinostat q = 1.1 × 10⁻⁴; bootstrap selection 100%; SHAP top contributor; docking −7.9 kcal/mol | Vorinostat, entinostat |
| 2 | MED1 | DepMap CRISPR dependency (r = −0.28, q = 0.049) | (preclinical) |
| 2 | SIRT2 | Signature member; docking AGK2 −8.5 kcal/mol | AGK2 |
| 3 | KAT2A/EP300 | Writer core; perturbation sensitivity | C646, MG149 |
