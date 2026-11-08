import re

s = open("build_gnd_manuscript.py", encoding="utf-8").read()

new_abs = '''# Abstract

**Background**: Transcript- and protein-level studies of the spliceosome yield contradictory prognostic signals in glioma. We asked whether glioblastoma (GBM) harbors a state marked by post-transcriptional spliceosome accumulation. **Methods**: We integrated single-cell (35,823 cells), spatial, and bulk transcriptomics with CPTAC proteomics, acetylomics, and phosphoproteomics, defining an expression-based Writer–Eraser Balance Score (WEBS). **Results**: A writer-dominant state showed proteome-wide spliceosome accumulation (107/111 proteins significantly up-regulated) that tracked WEBS at the protein (ρ = 0.79, 95% CI 0.68–0.87) but not the mRNA level (ρ = −0.03), although discordance magnitude was not spliceosome-specific under expression-matched permutation (p = 0.22). Spliceosome acetylation was enriched (parent-protein-corrected OR = 2.29), phosphosite signal changed in an abundance-dependent manner, and 4,387 alternative splicing events tracked WEBS, replicating in TCGA-LGG (83% sign-concordant). A six-gene signature was externally validated in three glioma cohorts (GLASS HR = 2.14, p = 2.0 × 10⁻⁶; C-index 0.62–0.67). Candidate vulnerabilities (MED1, SF3B1, HDAC2, SIRT2) were nominated computationally. **Conclusions**: A writer-dominant expression state is associated with post-transcriptional spliceosome accumulation, coordinated post-translational modification, and altered splicing output in GBM, nominating candidate vulnerabilities for experimental testing.

(196 words)'''

pattern = re.compile(r'new_abstract = """# Abstract.*?words\)"""', re.S)
s2 = pattern.sub(lambda m: 'new_abstract = """' + new_abs + '"""', s, count=1)
assert s2 != s, "abstract pattern not found"
open("build_gnd_manuscript.py", "w", encoding="utf-8").write(s2)
print("abstract patched to 196 words")
