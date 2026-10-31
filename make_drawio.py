# Generate Figure_1 and Figure_12 as draw.io (mxGraph) XML files
import os

OUT = r"E:\KIMI\GBM_WEBS_Project\figures_final"
os.makedirs(OUT, exist_ok=True)

ID = [0]
def nid():
    ID[0] += 1
    return f"c{ID[0]}"

def cell(id, value, x, y, w, h, style, edge=None, source=None, target=None, parent="1"):
    v = value.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    geom = f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/>'
    if edge:
        st = f'source="{source}" target="{target}"' if source and target else ""
        return (f'<mxCell id="{id}" value="{v}" style="{style}" '
                f'vertex="0" edge="1" parent="{parent}" {st}>{geom}</mxCell>')
    return f'<mxCell id="{id}" value="{v}" style="{style}" vertex="1" parent="{parent}">{geom}</mxCell>'

BOX = "rounded=1;whiteSpace=wrap;html=1;fillColor=#F5F5F5;strokeColor=#333333;fontColor=#1A1A1A;fontSize=12;fontFamily=Arial;"
def colorbox(hexcolor, font="#FFFFFF"):
    return (f"rounded=1;whiteSpace=wrap;html=1;fillColor={hexcolor};strokeColor=none;"
            f"fontColor={font};fontSize=12;fontFamily=Arial;fontStyle=1;")
ARROW = "edgeStyle=orthogonalEdgeStyle;rounded=0;html=1;endArrow=block;endFill=1;strokeColor=#333333;strokeWidth=2;fontFamily=Arial;"

LLS, LDTS, LHS = "#2E5C8A", "#9BB8D3", "#D95D5D"
PURPLE, TEAL, RED, ORANGE, GOLD = "#6B4C9A", "#3D9B8C", "#C0392B", "#D98E32", "#FFF8E1"

def doc(cells):
    return ('<mxfile host="app.diagrams.net"><diagram name="figure">'
            '<mxGraphModel dx="1100" dy="760" grid="1" gridSize="10" guides="1" '
            'tooltips="1" connect="1" arrows="1" fold="1" page="1" pageScale="1" '
            'pageWidth="1169" pageHeight="827" math="0" shadow="0">'
            '<root><mxCell id="0"/><mxCell id="1" parent="0"/>'
            + "".join(cells) + "</root></mxGraphModel></diagram></mxfile>")

# ============================== FIGURE 1 ==============================
c = []

# Title
c.append(cell(nid(), "A computational framework for discovering a therapeutically targetable malignant cell state in glioblastoma",
              80, 20, 900, 40,
              "text;html=1;align=center;verticalAlign=middle;fontSize=18;fontStyle=1;fontFamily=Arial;fontColor=#1A1A1A;"))

# Data sources (left)
c.append(cell(nid(), "DATA SOURCES", 60, 90, 220, 26,
              "text;html=1;align=left;verticalAlign=middle;fontSize=12;fontStyle=1;fontFamily=Arial;fontColor=#8A8A8A;"))
sources = [
    ("scRNA-seq (GSE138794)", "35,823 cells / 17 patients"),
    ("Spatial transcriptomics (GSE194329)", "6 Visium slices"),
    ("snATAC-seq (GSE138794)", "6 libraries"),
    ("Bulk cohorts (n = 8)", "TCGA · CGGA · GSE · GLASS"),
    ("CPTAC proteomics", "99 samples · 11,141 proteins"),
    ("DepMap / GDSC2", "CRISPR + drug screens"),
]
src_ids = []
for i, (t, sub) in enumerate(sources):
    y = 130 + i * 78
    sid = nid(); src_ids.append(sid)
    c.append(cell(sid, f"<b>{t}</b><br><font style='font-size:10px'>{sub}</font>", 40, y, 260, 62, BOX))

# WEBS center
c.append(cell(nid(), "<b>WEBS</b><br><font style='font-size:10px'>UCell(writers) − UCell(erasers)</font>", 430, 300, 180, 90,
              f"rounded=1;whiteSpace=wrap;html=1;fillColor=#FFFFFF;strokeColor={PURPLE};strokeWidth=2;fontColor={PURPLE};fontSize=16;fontFamily=Arial;fontStyle=1;"))
c.append(cell(nid(), "Writers<br>6 KATs", 430, 210, 85, 60, colorbox(PURPLE)))
c.append(cell(nid(), "Erasers<br>8 HDAC/SIRT", 525, 210, 85, 60, colorbox(TEAL)))
c.append(cell(nid(), "LLS", 430, 410, 55, 40, colorbox(LLS)))
c.append(cell(nid(), "LDTS", 492, 410, 55, 40, colorbox(LDTS)))
c.append(cell(nid(), "LHS", 554, 410, 55, 40, colorbox(LHS)))
webs_id = "c5"

# edges: sources -> webs
for i, sid in enumerate(src_ids):
    y = 130 + i * 78 + 31
    c.append(cell(nid(), "", 300, y, 130, 1, ARROW, edge=True, source=sid, target=webs_id))

# Tiers (right)
c.append(cell(nid(), "THREE-TIER ANALYSIS", 740, 90, 260, 26,
              "text;html=1;align=left;verticalAlign=middle;fontSize=12;fontStyle=1;fontFamily=Arial;fontColor=#8A8A8A;"))
tiers = [
    ("Discovery", "Cell-state landscape · regions<br>trajectory · hypoxia decoupling", LLS, 130),
    ("Mechanism", "CPTAC spliceosome activation<br>post-transcriptional · acetylation coupling", ORANGE, 300),
    ("Translation", "Six-gene signature · pharmacogenomics<br>docking · drug repositioning", RED, 470),
]
tier_ids = []
for label, body, col, y in tiers:
    tid = nid(); tier_ids.append(tid)
    c.append(cell(tid, f"<b>{label}</b><br><font style='font-size:10px;fontStyle:0'>{body}</font>",
                  720, y, 300, 110,
                  f"rounded=1;whiteSpace=wrap;html=1;fillColor=#FFFFFF;strokeColor={col};strokeWidth=2;fontColor={col};fontSize=13;fontFamily=Arial;"))
    c.append(cell(nid(), "", 610, 345, 110, 1, ARROW, edge=True, source=webs_id, target=tid))

# Key finding
kf = nid()
c.append(cell(kf, "<b>KEY FINDING</b><br>Post-transcriptional spliceosome activation<br><b><font color='#C0392B'>Protein ρ = +0.79  vs  mRNA ρ = −0.03</font></b>",
              720, 620, 300, 90,
              f"rounded=1;whiteSpace=wrap;html=1;fillColor={GOLD};strokeColor={RED};strokeWidth=2;fontColor=#1A1A1A;fontSize=11;fontFamily=Arial;"))

# Bottom bar
c.append(cell(nid(), "Reproducible pipeline · GitHub + Zenodo DOI · 30 scripts · WEBS standalone tool<br>CopyKAT · bootstrap · SHAP · parameter sensitivity · method benchmarking",
              40, 730, 700, 50,
              "text;html=1;align=left;verticalAlign=middle;fontSize=10;fontFamily=Arial;fontColor=#8A8A8A;"))

with open(os.path.join(OUT, "Figure_1.drawio"), "w", encoding="utf-8") as f:
    f.write(doc(c))

# ============================== FIGURE 12 ==============================
c = []
c.append(cell(nid(), "Working model of the epigenetic writer–eraser balance axis",
              80, 20, 900, 40,
              "text;html=1;align=center;verticalAlign=middle;fontSize=18;fontStyle=1;fontFamily=Arial;fontColor=#1A1A1A;"))

# Balance
c.append(cell(nid(), "Writers<br><font style='font-size:10px'>EP300 · KAT2A/B · KAT5 · KAT8 · CREBBP</font>", 60, 100, 180, 70,
              f"rounded=1;whiteSpace=wrap;html=1;fillColor=#FFFFFF;strokeColor={PURPLE};strokeWidth=2;fontColor={PURPLE};fontSize=13;fontFamily=Arial;fontStyle=1;"))
c.append(cell(nid(), "Erasers<br><font style='font-size:10px'>HDAC1-3 · SIRT1-3,6,7</font>", 280, 100, 180, 70,
              f"rounded=1;whiteSpace=wrap;html=1;fillColor=#FFFFFF;strokeColor={TEAL};strokeWidth=2;fontColor={TEAL};fontSize=13;fontFamily=Arial;fontStyle=1;"))
lhs_id = nid()
c.append(cell(lhs_id, "Writer-dominant state (LHS)<br><font style='font-size:10px'>IQR top 25% of WEBS</font>", 540, 100, 240, 70, colorbox(LHS)))
bal_w, bal_e = nid(), nid()
c.append(cell(nid(), "", 240, 135, 300, 1, ARROW, edge=True, target=lhs_id))

# Mechanism chain
c.append(cell(nid(), "POST-TRANSCRIPTIONAL MECHANISM", 60, 210, 400, 26,
              "text;html=1;align=left;verticalAlign=middle;fontSize=12;fontStyle=1;fontFamily=Arial;fontColor=#8A8A8A;"))
steps = [
    ("Non-histone<br>acetyltransferase<br>activity", 60),
    ("Splicing-factor<br>lysine acetylation<br><b>43 sites · OR = 7.36</b>", 300),
    ("Competitive block<br>of ubiquitination", 540),
    ("Protein<br>stabilization<br><b>ρ = 0.62</b>", 780),
]
step_ids = []
for label, x in steps:
    sid = nid(); step_ids.append(sid)
    c.append(cell(sid, label, x, 250, 200, 90, BOX))
for i in range(3):
    c.append(cell(nid(), "", 260, 295, 40, 1, ARROW, edge=True, source=step_ids[i], target=step_ids[i+1]))

# Core node
core = nid()
c.append(cell(core, "<b>SPLICEOSOME ACCUMULATION</b><br>107/111 proteins up-regulated · p = 3.4 × 10⁻³⁷ · protein ρ = +0.79 · mRNA ρ = −0.03",
              280, 390, 500, 80, colorbox(PURPLE)))
c.append(cell(nid(), "", 380, 295, 240, 95, ARROW, edge=True, source=step_ids[3], target=core))

# Downstream
c.append(cell(nid(), "CONSEQUENCES", 60, 500, 400, 26,
              "text;html=1;align=left;verticalAlign=middle;fontSize=12;fontStyle=1;fontFamily=Arial;fontColor=#8A8A8A;"))
cons = ["Enhanced<br>DNA repair", "Ribosome<br>biogenesis", "Survival &<br>proliferation", "Immunosuppressive<br>milieu (MIF–CD74<br>enriched in LLS)"]
cons_ids = []
for i, label in enumerate(cons):
    x = 60 + i * 240
    cid = nid(); cons_ids.append(cid)
    c.append(cell(cid, label, x, 540, 200, 80,
                  f"rounded=1;whiteSpace=wrap;html=1;fillColor=#FFFFFF;strokeColor=#8A8A8A;fontColor=#1A1A1A;fontSize=12;fontFamily=Arial;"))
    c.append(cell(nid(), "", 380, 470, 1, 70, ARROW, edge=True, source=core, target=cid))

# Therapeutics
c.append(cell(nid(), "THERAPEUTICS", 60, 650, 400, 26,
              "text;html=1;align=left;verticalAlign=middle;fontSize=12;fontStyle=1;fontFamily=Arial;fontColor=#8A8A8A;"))
c.append(cell(nid(), "<b>HDAC inhibitors</b><br>vorinostat · entinostat<br>HDAC2 docking −7.9 kcal/mol", 120, 690, 260, 90,
              f"rounded=1;whiteSpace=wrap;html=1;fillColor=#FFFFFF;strokeColor={RED};strokeWidth=2;fontColor={RED};fontSize=12;fontFamily=Arial;"))
c.append(cell(nid(), "<b>Spliceosome modulators</b><br>H3B-8800 · E7107<br>SF3B1 docking −6.7 kcal/mol", 620, 690, 260, 90,
              f"rounded=1;whiteSpace=wrap;html=1;fillColor=#FFFFFF;strokeColor={ORANGE};strokeWidth=2;fontColor={ORANGE};fontSize=12;fontFamily=Arial;"))
c.append(cell(nid(), "+", 480, 710, 40, 50,
              "text;html=1;align=center;verticalAlign=middle;fontSize=24;fontStyle=1;fontFamily=Arial;fontColor=#1A1A1A;"))

# Clinical anchor
c.append(cell(nid(), "Clinical anchor · six-gene risk score (HIF1A·SIRT2·MMP2·ITGA5·HDAC1·SPP1)<br>TCGA p = 0.0028 · CGGA-693 p = 3.4 × 10⁻¹⁹ · CGGA-325 p = 7.7 × 10⁻¹⁸",
              60, 810, 900, 50,
              "text;html=1;align=left;verticalAlign=middle;fontSize=10;fontFamily=Arial;fontColor=#8A8A8A;"))

with open(os.path.join(OUT, "Figure_12.drawio"), "w", encoding="utf-8") as f:
    f.write(doc(c))

print("drawio files created in", OUT)
