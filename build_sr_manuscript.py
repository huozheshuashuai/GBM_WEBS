import re

s = open("MANUSCRIPT_GND.md", encoding="utf-8").read()

# ---------- 1. unstructured single-paragraph abstract (~200 words) ----------
new_abs = """# Abstract

Transcript- and protein-level studies of the spliceosome yield contradictory prognostic signals in glioma. We asked whether glioblastoma (GBM) harbors a state marked by post-transcriptional spliceosome accumulation. We integrated single-cell (35,823 cells), spatial, and bulk transcriptomics with CPTAC proteomics, acetylomics, and phosphoproteomics, defining an expression-based Writer–Eraser Balance Score (WEBS). A writer-dominant state showed proteome-wide spliceosome accumulation (107/111 proteins significantly up-regulated) that tracked WEBS at the protein (ρ = 0.79, 95% CI 0.68–0.87) but not the mRNA level (ρ = −0.03), although discordance magnitude was not spliceosome-specific under expression-matched permutation (p = 0.22); per-protein mRNA-adjusted residuals confirmed systematic accumulation (103/109 positive, Fisher OR = 55.9). Spliceosome acetylation was enriched (parent-protein-corrected OR = 2.29), phosphosite signal changed in an abundance-dependent manner, and 4,387 alternative splicing events tracked WEBS, replicating in TCGA-LGG (83% sign-concordant). A six-gene signature was externally validated in three glioma cohorts (GLASS HR = 2.14, p = 2.0 × 10⁻⁶; C-index 0.62–0.67). Candidate vulnerabilities (MED1, SF3B1, HDAC2, SIRT2) were nominated computationally. A writer-dominant expression state is associated with post-transcriptional spliceosome accumulation, coordinated post-translational modification, and altered splicing output in GBM, nominating candidate vulnerabilities for experimental testing."""

m = re.search(r"# Abstract\n\n.*?\(196 words\)", s, re.S)
assert m, "abstract not found"
s = s[:m.start()] + new_abs + s[m.end():]

# ---------- 2. remove Highlights (SR does not use them) ----------
s = re.sub(r"# Highlights\n\n.*?(?=# Abbreviations)", "", s, flags = re.S)

# ---------- 3. reorder: Methods to the end (SR style) ----------
mm = re.search(r"(# 2\. Methods\n\n.*?)(# 3\. Results)", s, re.S)
results_end = re.search(r"# Acknowledgements", s)
methods_block = mm.group(1)
s = s[:mm.start()] + mm.group(2) + s[mm.end():]
# insert Methods after Conclusion, before Acknowledgements
s = s.replace("# Acknowledgements", methods_block + "\n# Acknowledgements")

# ---------- 4. Additional info per SR style ----------
addinfo = """# Acknowledgments

The authors thank the data-generating consortia (TCGA, CPTAC, CGGA, GLASS, DepMap, GDSC) and the original study authors.

# Author contributions

**Jia Chen**: Conceptualization, Data curation, Formal analysis, Investigation, Methodology, Software, Validation, Visualization, Writing – original draft. **Junhan Yang**: Validation, Writing – review & editing. **Bing Yang**: Writing – review & editing. **Tianyi Zhou**: Writing – review & editing. **Peijun Zhao**: Writing – review & editing.

# Competing interests

The authors declare no competing interests.

# Data availability

All datasets analyzed in this study are available from their respective repositories, subject to the access requirements and terms of use established by the original data providers: GSE138794, GSE194329, GSE131928, GSE7696, GSE83300 (GEO, open access); TCGA cohorts (GDC and the UCSC Xena GDC hub, STAR-TPM release); CGGA-693/325 (cgga.org.cn, registration and data-use agreement required); CPTAC-GBM (LinkedOmics); TCGASpliceSeq (bioinformatics.mdanderson.org); DepMap 24Q2 (depmap.org); GDSC2 (cancerrxgene.org). GLASS data were accessed via Synapse (syn17038081) under the consortium's data-use terms and are not redistributed. Analysis code and the WEBS implementation are publicly available at https://github.com/huozheshuashuai/GBM_WEBS; the submission version corresponds to the tagged release indicated in the repository README.

# Ethics declaration

This study involved only secondary analyses of publicly available, de-identified datasets. No new human participants, tissue samples, or identifiable data were collected. All datasets were accessed and analyzed in accordance with their respective data-use agreements.

# AI declaration

During the preparation of this work, the authors used generative AI tools (large language models) for language editing, code formatting, figure illustration assistance, and computational pipeline reconstruction assistance. All analyses were executed by the authors, and all code, outputs, and scientific conclusions were independently verified by the authors against the underlying data. The authors take full responsibility for the content of this publication.

# References"""

s = re.sub(r"# Ethics approval.*?# References", addinfo, s, flags = re.S)

open("MANUSCRIPT_SR.md", "w", encoding="utf-8").write(s)
print("MANUSCRIPT_SR.md written")
