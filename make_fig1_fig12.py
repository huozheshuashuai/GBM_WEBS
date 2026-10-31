# =============================================================================
# Publication-style Figure 1 (Study Design) & Figure 12 (Working Model)
# Vector-quality schematics via matplotlib
# =============================================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D
import os

OUT = r"E:\KIMI\GBM_WEBS_Project\figures_final"
os.makedirs(OUT, exist_ok=True)

# ---------- palette ----------
C = {
    "lls": "#2E5C8A", "ldts": "#7BA3C9", "lhs": "#E0665C",
    "writer": "#7B4EA3", "eraser": "#3D9B8C",
    "box": "#F5F5F5", "border": "#333333", "text": "#1A1A1A",
    "disc": "#2E5C8A", "mech": "#E08A3C", "trans": "#C0392B",
    "accent": "#B2182B", "gold": "#FFF8E1", "purple": "#6B4C9A",
}
plt.rcParams["font.family"] = "DejaVu Sans"

def rbox(ax, x, y, w, h, fc, ec=C["border"], lw=1.5, r=0.02, z=2):
    b = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.004,rounding_size={r}",
                       fc=fc, ec=ec, lw=lw, zorder=z)
    ax.add_patch(b)
    return b

def txt(ax, x, y, s, size=10, color=C["text"], weight="normal", ha="center", va="center", z=3):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight,
            ha=ha, va=va, zorder=z, linespacing=1.35)

def arrow(ax, x1, y1, x2, y2, color=C["border"], lw=2.5, style="-|>", ms=14, z=1):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=ms,
                        color=color, lw=lw, zorder=z)
    ax.add_patch(a)

# =============================================================================
# FIGURE 1 — Study Design
# =============================================================================
fig, ax = plt.subplots(figsize=(12, 8.6), dpi=300)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

txt(ax, 50, 97, "A computational framework for discovering a therapeutically targetable\nmalignant cell state in glioblastoma",
    size=14, weight="bold")

# ---- left: data sources ----
txt(ax, 15, 90, "Data sources", size=11.5, weight="bold")
sources = [
    ("scRNA-seq (GSE138794)", "35,823 cells / 17 patients", 81),
    ("Spatial (GSE194329)", "6 Visium slices", 72),
    ("snATAC-seq (GSE138794)", "6 libraries", 63),
    ("Bulk cohorts (n = 8)", "TCGA · CGGA · GSE · GLASS", 54),
    ("CPTAC proteomics", "99 samples · 11,141 proteins", 45),
    ("DepMap / GDSC2", "CRISPR + drug screens", 36),
]
for label, sub, y in sources:
    rbox(ax, 2, y - 3.2, 26, 7, C["box"])
    txt(ax, 15, y + 0.6, label, size=9.5, weight="bold")
    txt(ax, 15, y - 1.9, sub, size=8, color="#555555")

# ---- center: WEBS ----
rbox(ax, 36, 45, 28, 26, "white", ec=C["purple"], lw=2.5, r=0.04)
txt(ax, 50, 68, "WEBS", size=15, weight="bold", color=C["purple"])
rbox(ax, 38.5, 56.5, 10.5, 8, C["writer"], ec="none")
txt(ax, 43.75, 60.5, "Writers\n6 KATs", size=8, color="white", weight="bold")
rbox(ax, 51, 56.5, 10.5, 8, C["eraser"], ec="none")
txt(ax, 56.25, 60.5, "Erasers\n8 HDAC/SIRT", size=8, color="white", weight="bold")
txt(ax, 50, 52.5, "UCell(writer) − UCell(eraser)", size=8.5)
rbox(ax, 36.5, 46.5, 8.5, 5.5, C["lls"], ec="none"); txt(ax, 40.75, 49.2, "LLS", size=8.5, color="white", weight="bold")
rbox(ax, 45.75, 46.5, 8.5, 5.5, C["ldts"], ec="none"); txt(ax, 50, 49.2, "LDTS", size=8.5, color="white", weight="bold")
rbox(ax, 55, 46.5, 8.5, 5.5, C["lhs"], ec="none"); txt(ax, 59.25, 49.2, "LHS", size=8.5, color="white", weight="bold")
for _, _, y in sources:
    arrow(ax, 28, y - 0.1, 36, 58 if y > 58 else (52 if y > 46 else 47))

# ---- right: three tiers ----
tiers = [
    ("Discovery", "Cell-state landscape · Neftel states\nRegion & hypoxia decoupling · Trajectory", C["disc"], 78),
    ("Mechanism", "CPTAC spliceosome activation (107/111)\nPost-transcriptional accumulation\nAcetylation coupling (OR = 7.36) · MOFA2", C["mech"], 56),
    ("Translation", "Six-gene signature · 8-cohort survival\nDepMap MED1/SF3B1 · GDSC2 HDACi\nDrug repositioning · Docking", C["trans"], 32),
]
txt(ax, 81, 90, "Three-tier analysis", size=11.5, weight="bold")
for label, body, col, y in tiers:
    rbox(ax, 70, y - 8, 28, 15, col, ec="none", r=0.03)
    txt(ax, 84, y - 0.5, label, size=10.5, color="white", weight="bold")
    txt(ax, 84, y - 4.5, body, size=8, color="white")
arrow(ax, 64, 62, 70, 74); arrow(ax, 64, 58, 70, 52); arrow(ax, 64, 50, 70, 36)

# ---- key finding callout ----
rbox(ax, 62, 8, 36, 16, C["gold"], ec=C["accent"], lw=2, r=0.03)
txt(ax, 80, 21, "KEY FINDING", size=10, weight="bold", color=C["accent"])
txt(ax, 80, 15.5, "Post-transcriptional spliceosome activation", size=9, weight="bold")
txt(ax, 80, 11.5, "Protein ρ = +0.79  vs  mRNA ρ = −0.03", size=9.5, weight="bold", color=C["accent"])
arrow(ax, 80, 24, 72, 30, color=C["accent"], lw=2, style="-|>")

# ---- bottom: reproducibility ----
txt(ax, 34, 14, "Reproducible pipeline", size=10.5, weight="bold")
rbox(ax, 8, 4, 52, 7.5, C["box"])
txt(ax, 34, 7.8, "GitHub + Zenodo DOI · 30 analysis scripts · WEBS standalone tool\nCopyKAT · bootstrap · SHAP · parameter sensitivity · method benchmarking",
    size=8.5)

fig.savefig(os.path.join(OUT, "Figure_1_new.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "Figure_1_new.png"), bbox_inches="tight", dpi=300)
plt.close(fig)

# =============================================================================
# FIGURE 12 — Working Model
# =============================================================================
fig, ax = plt.subplots(figsize=(12, 9.5), dpi=300)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

txt(ax, 50, 97.5, "Working model: epigenetic writer–eraser balance drives post-transcriptional\nspliceosome accumulation in glioblastoma",
    size=13.5, weight="bold")

# ---- top: writer-eraser balance ----
rbox(ax, 8, 82, 18, 10, C["writer"], ec="none")
txt(ax, 17, 87, "Writers\nEP300 · KAT2A/B\nKAT5 · KAT8 · CREBBP", size=8.5, color="white", weight="bold")
rbox(ax, 30, 84, 18, 8, C["eraser"], ec="none")
txt(ax, 39, 88, "Erasers\nHDAC1-3 · SIRT1-3,6,7", size=8.5, color="white", weight="bold")
arrow(ax, 48, 87, 56, 87, lw=3)
rbox(ax, 56, 82, 36, 10, C["lhs"], ec="none")
txt(ax, 74, 87, "WRITER-DOMINANT STATE (LHS)\nIQR top 25% of WEBS", size=10, color="white", weight="bold")

# ---- mechanism chain ----
txt(ax, 50, 76, "Post-transcriptional mechanism", size=11, weight="bold")
chain = [
    ("Non-histone\nacetyltransferase\nactivity", "#F2C14E", 8),
    ("Splicing-factor\nlysine acetylation\n43 sites · OR = 7.36", "#F5E1A4", 31),
    ("Blocks lysine\nubiquitination\n(competition)", "#BFE3C0", 54),
    ("Protein\nstabilization\nρ = 0.62", "#E08A3C", 77),
]
for label, col, x in chain:
    rbox(ax, x, 60, 16, 12, col, ec=C["border"])
    txt(ax, x + 8, 66, label, size=8, weight="bold")
for i in range(3):
    arrow(ax, chain[i][2] + 16, 66, chain[i + 1][2], 66, lw=2.5)

# ---- core node ----
rbox(ax, 20, 42, 60, 12, C["purple"], ec="none", r=0.03)
txt(ax, 50, 50, "SPLICEOSOME ACCUMULATION", size=13, color="white", weight="bold")
txt(ax, 50, 45.5, "107/111 proteins up-regulated · enrichment p = 3.4 × 10⁻³⁷\nprotein ρ = +0.79 · mRNA ρ = −0.03 (post-transcriptional)",
    size=8.5, color="white")
arrow(ax, 50, 60, 50, 54, lw=3)

# ---- downstream ----
down = [
    ("Enhanced\nDNA repair", "#6C8EBF", 10),
    ("Ribosome\nbiogenesis", "#9673A6", 32),
    ("Survival &\nproliferation", "#D79B9B", 54),
    ("Immunosuppressive\nmilieu (MIF–CD74\nenriched in LLS)", "#E8A87C", 76),
]
txt(ax, 50, 38, "Downstream consequences", size=11, weight="bold")
for label, col, x in down:
    rbox(ax, x, 22, 18, 11, col, ec="none")
    txt(ax, x + 9, 27.5, label, size=8, color="white", weight="bold")
    arrow(ax, 50, 42, x + 9, 33.5, lw=1.8, color="#777777")

# ---- therapeutics ----
txt(ax, 50, 16, "Therapeutic strategies", size=11, weight="bold")
rbox(ax, 12, 4, 32, 9, "#C0392B", ec="none")
txt(ax, 28, 8.5, "HDAC inhibitors\nVorinostat · Entinostat\nHDAC2 docking −7.9 kcal/mol", size=8, color="white", weight="bold")
rbox(ax, 56, 4, 32, 9, "#D35400", ec="none")
txt(ax, 72, 8.5, "Spliceosome modulators\nH3B-8800 · E7107\nSF3B1 docking −6.7 kcal/mol", size=8, color="white", weight="bold")
arrow(ax, 30, 22, 28, 13, lw=2, color="#777777")
arrow(ax, 70, 22, 72, 13, lw=2, color="#777777")
txt(ax, 50, 8.5, "+", size=22, weight="bold", color=C["border"])

# ---- clinical anchor ----
rbox(ax, 2, -0.5, 96, 4.5, C["box"])
txt(ax, 50, 2.2, "Clinical anchor: six-gene risk score (HIF1A·SIRT2·MMP2·ITGA5·HDAC1·SPP1)  |  TCGA p=0.0028 · CGGA-693 p=3.4e-19 · CGGA-325 p=7.7e-18",
    size=8, weight="bold")

ax.set_ylim(-1, 100)
fig.savefig(os.path.join(OUT, "Figure_12_new.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "Figure_12_new.png"), bbox_inches="tight", dpi=300)
plt.close(fig)

print("Figure 1 & 12 generated ->", OUT)
