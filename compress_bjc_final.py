import re

s = open("MANUSCRIPT_BJC.md", encoding="utf-8").read()

# ---- 3.9 compression ----
old39 = re.search(r"## 3\.9 A bootstrap-stable six-gene signature.*?## 3\.10", s, re.S)
new39 = """## 3.9 A bootstrap-stable six-gene signature, externally validated in three independent cohorts including GLASS

The six-gene signature (HIF1A, SIRT2, MMP2, ITGA5, HDAC1, SPP1; fixed coefficients) stratified survival in TCGA-GBM (log-rank p = 0.0028), CGGA-693 (p = 3.4 × 10⁻¹⁹), CGGA-325 (p = 7.7 × 10⁻¹⁸), and GLASS (HR = 2.14, p = 2.0 × 10⁻⁶; Fig. 6). Two additional cohorts showed concordant directions without significance (GSE7696, GSE83300; Fig. S11); attenuation in IDH-wild-type subsets is reported explicitly as a boundary (Fig. S28). Exploiting the longitudinal design of GLASS, higher primary-tumor WEBS was associated with a longer observed interval to first recurrence (ρ = 0.33, p = 1.1 × 10⁻⁶, 213 recurrence-interval records from 173 patients; confirmed per-patient, ρ = 0.31; Fig. 6e); this exploratory association does not account for censoring or treatment variables and is subject to re-sampling selection bias. WEBS itself did not change systematically at recurrence (179 pairs from 148 patients; Fig. S14).

The signature's full performance profile (Table S_signature_performance: per-SD HRs 1.14–1.73; C-indices 0.62–0.67; 3-year AUCs 0.67–0.76) indicates a moderate, auxiliary stratification marker. Adjustment for clinical covariates retained significance in CGGA-325 (p = 0.0011) but attenuated in CGGA-693 (p = 0.182) — reported transparently. Exploratory risk × treatment interactions (uncorrected) are shown in Fig. S12, and joint stratification with MGMT methylation improved discrimination (p < 10⁻¹³; Fig. S13). Given substantial between-cohort heterogeneity (I² = 91.4%), the random-effects meta-analysis was the primary pooled estimate (HR = 1.35, p = 0.0028; common-effect HR = 1.23 as sensitivity; Fig. S15). Proportional-hazards tests (cox.zph) were acceptable in TCGA-GBM and GLASS (global p = 0.68, 0.20) but violated in the two CGGA cohorts (p = 0.020, 0.0036), noted as a limitation. Bootstrap (1,000×) confirmed feature-selection stability (HDAC1 in 99.9%; Fig. S5); SHAP identified HDAC1 as the dominant protective contributor (Fig. S6), and ablation analysis showed a distributed rather than single-gene-driven signature (Table S_ablation). An RPPA lactate-metabolism score stratified survival at the protein level (p = 0.016, n = 71).

## 3.10"""
if old39:
    s = s[:old39.start()] + new39 + s[old39.end():]

# ---- 2.5 snATAC compression ----
old25 = re.search(r"## 2\.5 snATAC-seq\n\n.*?\n\n## 2\.6", s, re.S)
new25 = """## 2.5 snATAC-seq

snATAC-seq libraries from the same study (GSE138794; six libraries: GSM4119513–GSM4119518) were analyzed. Per GEO annotation, these comprise two IDH-wild-type GBM, two IDH-mutant astrocytomas, and two oligodendrogliomas (one recurrent) — an exploratory mixed-histology glioma analysis, not GBM-specific. Fragment files were not deposited; analysis was performed at peak-matrix level (Signac [20]). Peaks are GRCh38 (verified by TSS-overlap test); promoter accessibility was aggregated over TSS ± 2 kb of writer/eraser genes (EnsDb.Hsapiens.v86) and correlated with patient-median scRNA WEBS.
"""
if old25:
    s = s[:old25.start()] + new25 + s[old25.end():]

# ---- 2.10 compression ----
old210 = re.search(r"## 2\.10 Prognostic modeling and longitudinal analysis\n\n.*?\n\n## 2\.11", s, re.S)
new210 = """## 2.10 Prognostic modeling and longitudinal analysis

LASSO-Cox [21] (glmnet, 10-fold CV, lambda.min) on 29 pre-specified lactylation/hypoxia/invasion candidates was trained in TCGA-GBM (LinkedOmics, n = 144), yielding the six-gene signature (HIF1A, SIRT2, MMP2, ITGA5, HDAC1, SPP1) with coefficients locked thereafter [22]; the candidate pool was defined a priori (not re-selected within folds), and nested cross-validation was not performed (noted as a limitation). Performance was evaluated in three independent external cohorts (CGGA-693, CGGA-325, GLASS); GSE7696 and GSE83300 were additional underpowered cohorts, and IDH-wild-type analyses were subgroup analyses only. The TCGA-GDC cohort (containing all 144 training patients) was used only for boundary assessment. Risk scores were computed from log2-normalized expression with fixed coefficients in every cohort (no re-standardization or refitting) and dichotomized at the cohort-specific median. Bootstrap (1,000× full-pipeline) selection frequencies; SHAP (kernelshap) on the Cox linear predictor; RPPA lactate-metabolism score (n = 71). GLASS longitudinal analysis: 179 primary–recurrence pairs (148 patients) by paired tests; 213 recurrence-interval records (173 patients) with a first-record-per-patient sensitivity analysis.
"""
if old210:
    s = s[:old210.start()] + new210 + s[old210.end():]

# ---- 3.6 compression ----
old36 = re.search(r"## 3\.6 The accumulated spliceosome shows signatures of functional engagement\n\n.*?\n\n## 3\.7", s, re.S)
new36 = """## 3.6 The accumulated spliceosome shows signatures of functional engagement

Accumulation of an enzymatic machine does not itself imply altered activity, so we tested functional output directly. Alternative splicing analysis (TCGASpliceSeq, TCGA-GBM) identified 4,387 of 28,071 events (15.6%) significantly associated with bulk WEBS (q < 0.05, |ρ| > 0.2; Fig. 5a,b). Genes harboring these events were enriched for RNA splicing and mRNA processing (adjusted p = 1.2 × 10⁻³; Fig. 5c). The association was not fully explained by immune-content confounding (91% of events retained after immune-marker partial correlation), and 97% (4,243/4,387) retained significance after simultaneous adjustment for immune, proliferation, mesenchymal, proneural, and hypoxia programs — indicating the association is not driven by cell composition or molecular-subtype programs. The events replicated cross-cohort in TCGA-LGG (83% sign-concordant of 4,285 testable; effect ρ = 0.70; Fig. 5d).

Independently, the CPTAC phosphoproteome showed more spliceosome phosphosite signal in LHS (198/319 sites, 197 increased; enrichment OR = 3.26). After parent-protein correction, however, most remaining significant sites decreased (118/124), including SR-protein RS domains (27/28 SRSF sites; Fig. 4b) — the raw signal is largely abundance-driven while the abundance-adjusted signal is reduced. Splicing-kinase abundance tracked WEBS (CLK3, CDK9/12/13; Fig. 4c). Thus the accumulated machinery is not simply "hyper-phosphorylated"; its phosphorylation stoichiometry is altered in a site- and layer-dependent manner.

Finally, connectivity-based sensitivity analysis (GENIE3, unsigned weights) showed spliceosome genes to be more sensitive than other genes for all 14 writer/eraser regulators (ratios 1.31–1.53; regulator-level p = 6.1 × 10⁻⁵; Fig. 4d), indicating disproportionate wiring to acetylation-regulatory input. Together, these three layers — altered splicing output, abundance-dependent phosphorylation changes, and network sensitivity — support a functionally engaged spliceosome program, although direct perturbation experiments are required for causal claims.

## 3.7"""
if old36:
    s = s[:old36.start()] + new36 + s[old36.end():]

# ---- 4.4 compression ----
old44 = re.search(r"## 4\.4 Therapeutic prioritization\n\n.*?\n\n## 4\.5", s, re.S)
new44 = """## 4.4 Therapeutic prioritization

Table 2 ranks candidate targets by convergent evidence, where Tier 1 denotes priority for experimental validation rather than validated efficacy. Two candidate mechanisms emerge — spliceosome inhibition downstream (H3B-8800, E7107) and HDAC inhibition at the balance level (vorinostat, entinostat). Network-based sensitivity analysis and docking provided complementary computational hypotheses for HDAC2 and SIRT2, whereas GDSC2 independently identified associations between WEBS and sensitivity to selected HDAC inhibitors. These findings prioritize candidates for experimental testing but do not establish target-specific dependency or drug efficacy. The GLASS recurrence-interval association is hypothesis-generating and subject to longitudinal-sampling selection bias.

## 4.5"""
if old44:
    s = s[:old44.start()] + new44 + s[old44.end():]

# ---- 4.5 compression ----
old45 = re.search(r"## 4\.5 Experimental validation roadmap\n\n.*?\n\n## 4\.6", s, re.S)
new45 = """## 4.5 Experimental validation roadmap

1. **Targeted proteomics**: SF3B1/PRPF8 protein vs mRNA after EP300 overexpression or HDAC1/2 knockdown.
2. **Acetylation rescue**: mutate top-coupled lysines (PRPF40B, HNRNPA2B1) and measure half-life ± vorinostat.
3. **Splicing-output validation**: RT-PCR of top WEBS-associated exons (IFI44, DCTN3, MYL6) after HDAC/CLK inhibition.
4. **Functional validation**: H3B-8800 and vorinostat dose–response in high- vs low-WEBS GBM lines or organoids.

## 4.6"""
if old45:
    s = s[:old45.start()] + new45 + s[old45.end():]

open("MANUSCRIPT_BJC.md", "w", encoding="utf-8").write(s)

parts = re.split(r"(?m)^(# .+)$", s)
cur = ""; counts = {}
for p in parts:
    if p.startswith("# "):
        cur = p.strip(); counts.setdefault(cur, 0)
    else:
        counts[cur] = counts.get(cur, 0) + len(re.findall(r"[\w’'-]+", p))
body = sum(counts.get(k, 0) for k in ["# 1. Introduction", "# 2. Methods", "# 3. Results", "# 4. Discussion", "# 5. Conclusion"])
print("body:", body)
