import re

src = open("MANUSCRIPT_v2.md", encoding="utf-8").read()

# ---------- 1. Abstract: structured, <=250 words ----------
new_abstract = """# Abstract

**Background**: Spliceosome dysregulation is a cancer hallmark, but reports on splicing machinery and prognosis in glioma are contradictory across measurement layers. We asked whether glioblastoma (GBM) harbors a malignant cell state marked by post-transcriptional spliceosome accumulation.

**Methods**: We integrated single-cell (35,823 cells), spatial, and bulk transcriptomics with CPTAC proteomics, acetylomics, and phosphoproteomics, defining a Writer–Eraser Balance Score (WEBS). Functional output was assessed by transcriptome-wide alternative splicing analysis (TCGASpliceSeq) and in-silico network perturbation. Clinical associations were tested in TCGA-GBM, CGGA-693, CGGA-325, GLASS, and 33 pan-cancer cohorts.

**Results**: A writer-dominant state showed proteome-wide spliceosome accumulation (107/111 proteins up-regulated, robust across stratification thresholds) that tracked WEBS at the protein (Spearman ρ = 0.79, 95% CI 0.68–0.87) but not the mRNA level (ρ = −0.03; bootstrap Δρ 95% CI 0.60–1.03). Spliceosome proteins showed enriched acetylation (OR = 7.36; parent-protein-corrected OR = 2.29) and coordinated hyper-phosphorylation (198/319 sites), with elevated splicing kinases (CLK3, CDK9/12/13). WEBS associated with 4,387 alternative splicing events that replicated in TCGA-LGG (83% sign-concordant, effect ρ = 0.70). A six-gene signature derived in TCGA-GBM was externally validated in CGGA-693 (p = 3.4 × 10⁻¹⁹), CGGA-325 (p = 7.7 × 10⁻¹⁸), and GLASS (HR = 2.14, p = 2.0 × 10⁻⁶); primary-tumor WEBS was associated with a longer observed recurrence interval. Computational evidence nominated candidate vulnerabilities (MED1, SF3B1 dependency; HDAC-inhibitor sensitivity; vorinostat–HDAC2 docking −7.9 kcal/mol).

**Conclusions**: A writer-dominant GBM state is associated with broad, post-transcriptional spliceosome accumulation accompanied by coordinated post-translational modification and altered splicing output, nominating candidate therapeutic vulnerabilities for experimental testing.

(249 words)"""

src = re.sub(r"# Abstract\n\n.*?\n\n\(\d+ words\)", new_abstract, src, flags = re.S)

# ---------- 2. Figure reference remapping (longest first) ----------
# old main -> new mapping
# Fig. 13d -> Fig. 4b ; 13e -> 4c ; 13f -> 4d ; 13a,b -> 5a ; 13c -> 5c
# Fig. 14a -> Fig. 6d ; 14b -> Fig. 7a ; 14c,d -> Fig. S32
# Fig. 10 -> Fig. 6 (KM) ; Fig. 11 -> Fig. 7c,d ; Fig. 12 -> Fig. S33
# Fig. 9 -> S31 ; Fig. 8 -> S30 ; Fig. 7 -> S29 ; Fig. 6 -> S28 ; Fig. 5 -> S27 ; Fig. 4 -> S26 ; Fig. 3 -> S25
reps = [
    ("Fig. 13a,b", "Fig. 5a"), ("Fig. 13d", "Fig. 4b"), ("Fig. 13e", "Fig. 4c"),
    ("Fig. 13f", "Fig. 4d"), ("Fig. 13c", "Fig. 5c"),
    ("Fig. 14a", "Fig. 6d"), ("Fig. 14b", "Fig. 7a"), ("Fig. 14c", "Fig. S32"),
    ("Fig. 14d", "Fig. S32"),
    ("Fig. 10)", "Fig. 6)"), ("(Fig. 10)", "(Fig. 6)"),
    ("Fig. 11)", "Fig. 7c,d)"),
    ("(Fig. 12)", "(Fig. S33)"),
    ("Fig. 9)", "Fig. S31)"), ("Fig. 8,", "Fig. S30,"),
    ("Fig. 7,", "Fig. S29,"),
    ("Fig. 6),", "Fig. S28),"), ("(Fig. 6)", "(Fig. S28)"),
    ("Fig. 5)", "Fig. S27)"),
    ("Fig. 4a", "Fig. S26a"), ("Fig. 4b", "S26b"), ("Fig. 4c,d", "Fig. S26c,d"),
    ("Fig. 2, Fig. 3", "Fig. 2, Fig. S25"),
    ("(all p > 0.4; Fig. 8)", "(all p > 0.4; Fig. S30)"),
    ("is summarized in Fig. 9.", "is summarized in Fig. S31."),
    ("(Fig. 2, Fig. S25)", "(Fig. 2, Fig. S25; external-cohort panels in Fig. S24)"),
]
for old, new in reps:
    src = src.replace(old, new)

# ---------- 3. Figure legends: rewrite for 7F ----------
new_legends = """# Figure Legends

**Figure 1. Study design and analytical workflow.** Data sources, WEBS computation, and the analysis tiers (discovery, functional validation, clinical association, candidate vulnerability nomination).

**Figure 2. The Writer–Eraser Balance Score (WEBS) landscape.** (a) Cell-type annotation of 35,823 GBM cells. (b) WEBS projected onto the UMAP. (c) WEBS by IQR group (LLS/LDTS/LHS; Kruskal–Wallis p < 10⁻³⁰⁰).

**Figure 3. Protein-level spliceosome accumulation is associated with the writer-dominant state (CPTAC-GBM).** (a) Proteome differential abundance (high- vs low-WEBS); 107/111 spliceosome proteins up-regulated (purple). (b) Protein-level WEBS vs spliceosome protein score (ρ = 0.79, 95% CI 0.68–0.87). (c) Per-gene RNA- vs protein-layer correlation with WEBS (107/109 above diagonal). (d) Robustness across stratification thresholds (109–111/111 up-regulated).

**Figure 4. Coordinated post-translational modification of the spliceosome.** (a) Top differential spliceosome acetylation sites, samples ordered by WEBS (enrichment OR = 7.36; abundance-corrected OR = 2.29). (b) Phosphosite volcano (LHS vs LLS); 198/319 spliceosome sites significant (corrected 124/319). (c) Splicing-kinase abundance vs WEBS. (d) GENIE3-based connectivity sensitivity analysis: relative absolute predicted shift of spliceosome vs other genes per regulator perturbation (dashed line = ratio 1; purple = writers, teal = erasers).

**Figure 5. Altered alternative splicing output.** (a) 28,071 AS events vs bulk WEBS (TCGA-GBM); 4,387 significant. (b) Event-type composition. (c) GO enrichment (top: RNA splicing). (d) Cross-cohort replication in TCGA-LGG (83% sign-concordant; effect ρ = 0.70).

**Figure 6. Clinical validation of the six-gene signature.** Kaplan–Meier curves with number-at-risk tables: (a) TCGA-GBM (training), (b) CGGA-693, (c) CGGA-325, (d) GLASS primary tumors. In-panel p values are log-rank; subgroup HRs with 95% CI are from Cox models. Risk groups were defined by the cohort-specific median.

**Figure 7. Longitudinal association and candidate vulnerabilities.** (a) Primary-tumor WEBS vs observed interval to first recurrence (GLASS, 213 pairs/173 patients; ρ = 0.33; exploratory). (b) WEBS in paired primary–recurrence samples. (c) Top WEBS-sensitizing compounds (GDSC2, q < 0.1). (d) DepMap dependency correlations for MED1, SF3B1, SRSF1, HDAC2, SIRT2.

# Tables"""

src = re.sub(r"# Figure Legends\n\n.*?# Tables", new_legends, src, flags = re.S)

# ---------- 4. Tables section: keep only Table 1; Table 2 -> supplementary ----------
src = src.replace("**Table 2. Computationally prioritized candidate therapeutic vulnerabilities.**",
                  "**(Table 2 is provided in the Supplementary Materials.)**")
# remove Table 2 markdown block from main tables section
src = re.sub(r"\n\| Tier \| Candidate \| Key evidence \| Tool compounds \|.*?(?=\n\n# Figures|\Z)",
             "\n", src, flags = re.S)

# ---------- 5. Figures section: new paths ----------
fig_section = """# Figures

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```
![Figure 1](submission_GND/Figures/jpg/Figure_1.jpg)

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```
![Figure 2](submission_GND/Figures/jpg/Figure_2.jpg)

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```
![Figure 3](submission_GND/Figures/jpg/Figure_3.jpg)

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```
![Figure 4](submission_GND/Figures/jpg/Figure_4.jpg)

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```
![Figure 5](submission_GND/Figures/jpg/Figure_5.jpg)

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```
![Figure 6](submission_GND/Figures/jpg/Figure_6.jpg)

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```
![Figure 7](submission_GND/Figures/jpg/Figure_7.jpg)
"""
src = re.sub(r"# Figures\n\n.*", fig_section, src, flags = re.S)

# ---------- 6. Results section 3.10 title fix (pharmaco figure ref) ----------
src = src.replace("(Fig. 11)", "(Fig. 7c,d)")
src = src.replace("(Fig. 7c,d). Notably", "(Fig. 7c,d). Notably")

open("MANUSCRIPT_GND.md", "w", encoding="utf-8").write(src)
print("MANUSCRIPT_GND.md written,", len(src), "chars")
