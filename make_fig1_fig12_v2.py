# =============================================================================
# Figure 1 & 12 — minimal editorial style (clean, airy, journal-grade schematic)
# =============================================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib import font_manager
import os

OUT = r"E:\KIMI\GBM_WEBS_Project\figures_final"
os.makedirs(OUT, exist_ok=True)

for fp in [r"C:\Windows\Fonts\arial.ttf", r"C:\Windows\Fonts\arialbd.ttf"]:
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.family"] = "Arial"

INK   = "#222222"   # near-black text
GREY  = "#8A8A8A"   # secondary text
LINE  = "#B9B9B9"   # hairline
LLS, LDTS, LHS = "#2E5C8A", "#9BB8D3", "#D95D5D"
PURPLE, TEAL, RED, ORANGE = "#6B4C9A", "#3D9B8C", "#C0392B", "#D98E32"

def rbox(ax, x, y, w, h, fc="white", ec=INK, lw=1.0, r=0.012):
    b = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.002,rounding_size={r}",
                       fc=fc, ec=ec, lw=lw, zorder=2)
    ax.add_patch(b)

def T(ax, x, y, s, size=10.5, color=INK, weight="normal", ha="center", va="center"):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight,
            ha=ha, va=va, zorder=3, linespacing=1.45)

def arr(ax, x1, y1, x2, y2, color=INK, lw=1.4):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=13, color=color, lw=lw, zorder=1))

# =============================================================================
# FIGURE 1
# =============================================================================
fig, ax = plt.subplots(figsize=(11.2, 7.2), dpi=300)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

T(ax, 50, 96, "Study design and analytical workflow", size=15, weight="bold")

# left column: data sources
T(ax, 16, 89, "DATA SOURCES", size=9.5, color=GREY, weight="bold")
sources = [
    "scRNA-seq · 35,823 cells",
    "Spatial transcriptomics · 6 slices",
    "snATAC-seq · 6 libraries",
    "Bulk cohorts · 8 cohorts",
    "CPTAC proteomics · 99 samples",
    "DepMap · GDSC2 screens",
]
for i, s in enumerate(sources):
    y = 82 - i * 8.2
    rbox(ax, 3, y - 3.0, 26, 6.2)
    T(ax, 16, y, s, size=9.5)

# center: WEBS
rbox(ax, 39, 46, 22, 18, fc="white", ec=PURPLE, lw=1.6, r=0.02)
T(ax, 50, 60.5, "WEBS", size=13, weight="bold", color=PURPLE)
T(ax, 50, 55.5, "UCell(writers) − UCell(erasers)\n6 KATs vs 8 HDAC/SIRT", size=9)
for j, (lab, col) in enumerate([("LLS", LLS), ("LDTS", LDTS), ("LHS", LHS)]):
    rbox(ax, 40.5 + j * 6.7, 47.5, 6, 4.5, fc=col, ec="none")
    T(ax, 43.5 + j * 6.7, 49.7, lab, size=8.5, color="white", weight="bold")
T(ax, 50, 41.5, "LLS · LDTS · LHS  (IQR stratification)", size=8, color=GREY)

for i in range(len(sources)):
    y = 82 - i * 8.2
    arr(ax, 29, y, 39, 56 if y > 64 else (55 if y > 52 else 50))

# right column: tiers
T(ax, 82, 89, "THREE-TIER ANALYSIS", size=9.5, color=GREY, weight="bold")
tiers = [
    ("Discovery", "Cell-state landscape · regions\ntrajectory · hypoxia decoupling", LLS, 76),
    ("Mechanism", "CPTAC spliceosome activation\npost-transcriptional · acetylation", ORANGE, 55),
    ("Translation", "six-gene signature · pharmacogenomics\ndocking · drug repositioning", RED, 34),
]
for label, body, col, y in tiers:
    rbox(ax, 70, y - 6.5, 27, 13, fc="white", ec=col, lw=1.6)
    T(ax, 83.5, y + 3.2, label, size=10.5, weight="bold", color=col)
    T(ax, 83.5, y - 1.5, body, size=8.5)
arr(ax, 61, 58, 70, 70); arr(ax, 61, 55, 70, 55); arr(ax, 61, 52, 70, 40)

# key finding strip
rbox(ax, 64, 12, 34, 14, fc="white", ec=RED, lw=1.8)
T(ax, 81, 22.5, "KEY FINDING", size=8.5, weight="bold", color=RED)
T(ax, 81, 17.2, "Post-transcriptional spliceosome activation\nprotein ρ = +0.79 · mRNA ρ = −0.03", size=9)

# bottom bar
ax.plot([3, 97], [7.5, 7.5], color=LINE, lw=0.8, zorder=1)
T(ax, 50, 3.5, "Reproducible pipeline · GitHub + Zenodo DOI · 30 scripts · WEBS standalone tool\nCopyKAT · bootstrap · SHAP · parameter sensitivity · method benchmarking",
    size=8.5, color=GREY)

fig.savefig(os.path.join(OUT, "Figure_1_new.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "Figure_1_new.png"), bbox_inches="tight", dpi=300)
plt.close(fig)

# =============================================================================
# FIGURE 12
# =============================================================================
fig, ax = plt.subplots(figsize=(11.2, 8.2), dpi=300)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

T(ax, 50, 96.5, "Working model of the writer–eraser balance axis", size=15, weight="bold")

# top: balance
T(ax, 25, 90, "REGULATION", size=9, color=GREY, weight="bold")
rbox(ax, 8, 80, 15, 8, fc="white", ec=PURPLE, lw=1.5)
T(ax, 15.5, 84, "Writers\n6 KATs", size=9, weight="bold", color=PURPLE)
rbox(ax, 27, 80, 15, 8, fc="white", ec=TEAL, lw=1.5)
T(ax, 34.5, 84, "Erasers\n8 HDAC/SIRT", size=9, weight="bold", color=TEAL)
arr(ax, 42, 84, 50, 84)
rbox(ax, 50, 80, 20, 8, fc=LHS, ec="none")
T(ax, 60, 84, "Writer-dominant state (LHS)", size=9.5, color="white", weight="bold")

# mechanism chain
T(ax, 50, 73.5, "POST-TRANSCRIPTIONAL MECHANISM", size=9, color=GREY, weight="bold")
steps = [
    "non-histone\nacetyltransferase\nactivity",
    "splicing-factor\nlysine acetylation\n(43 sites, OR = 7.36)",
    "competitive\nblock of\nubiquitination",
    "protein\nstabilization\n(ρ = 0.62)",
]
for i, s in enumerate(steps):
    x = 7 + i * 23.5
    rbox(ax, x, 59, 19, 10.5)
    T(ax, x + 9.5, 64.2, s, size=8.5)
    if i < 3:
        arr(ax, x + 19, 64, x + 23.5, 64)

# core
rbox(ax, 24, 42, 52, 11, fc=PURPLE, ec="none", r=0.02)
T(ax, 50, 50.2, "SPLICEOSOME ACCUMULATION", size=12, color="white", weight="bold")
T(ax, 50, 45.6, "107/111 proteins up · p = 3.4 × 10⁻³⁷ · protein ρ = +0.79 · mRNA ρ = −0.03",
    size=8.5, color="white")
arr(ax, 50, 59, 50, 53.5)

# downstream
T(ax, 50, 37.5, "CONSEQUENCES", size=9, color=GREY, weight="bold")
cons = ["enhanced\nDNA repair", "ribosome\nbiogenesis", "survival &\nproliferation", "immunosuppressive\nmilieu (MIF–CD74\nin LLS)"]
for i, s in enumerate(cons):
    x = 7 + i * 23.5
    rbox(ax, x, 25, 19, 9, fc="white", ec=GREY, lw=1.0)
    T(ax, x + 9.5, 29.5, s, size=8.5)
    arr(ax, 50, 42, x + 9.5, 34.5, color=GREY, lw=1.0)

# therapeutics
T(ax, 50, 19.5, "THERAPEUTICS", size=9, color=GREY, weight="bold")
rbox(ax, 14, 7.5, 30, 9.5, fc="white", ec=RED, lw=1.6)
T(ax, 29, 13.5, "HDAC inhibitors", size=10, weight="bold", color=RED)
T(ax, 29, 10, "vorinostat · entinostat\nHDAC2 docking −7.9 kcal/mol", size=8.5)
rbox(ax, 56, 7.5, 30, 9.5, fc="white", ec=ORANGE, lw=1.6)
T(ax, 71, 13.5, "Spliceosome modulators", size=10, weight="bold", color=ORANGE)
T(ax, 71, 10, "H3B-8800 · E7107\nSF3B1 docking −6.7 kcal/mol", size=8.5)
T(ax, 50, 12, "+", size=20, weight="bold")
arr(ax, 29, 25, 29, 17.5, color=GREY, lw=1.0)
arr(ax, 71, 25, 71, 17.5, color=GREY, lw=1.0)

# clinical anchor
ax.plot([3, 97], [4.5, 4.5], color=LINE, lw=0.8, zorder=1)
T(ax, 50, 1.2, "Clinical anchor · six-gene risk score (HIF1A·SIRT2·MMP2·ITGA5·HDAC1·SPP1)\nTCGA p = 0.0028 · CGGA-693 p = 3.4 × 10⁻¹⁹ · CGGA-325 p = 7.7 × 10⁻¹⁸",
    size=8.5, color=GREY)

fig.savefig(os.path.join(OUT, "Figure_12_new.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "Figure_12_new.png"), bbox_inches="tight", dpi=300)
plt.close(fig)

print("v2 done")
