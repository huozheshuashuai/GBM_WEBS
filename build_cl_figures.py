# Compose 6 main figures for Cancer Letters from existing analysis PDFs
import fitz, os, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

BASE = r"E:/KIMI/GBM_WEBS_Project"
OUT = os.path.join(BASE, "submission_CancerLetters", "02_Figures")
os.makedirs(OUT, exist_ok=True)

def compose(panels, out_pdf, cols, cell_w=560, cell_h=430, title=None):
    """panels: list of (pdf_path, label). Grid layout."""
    rows = (len(panels) + cols - 1) // cols
    W, H = cell_w * cols, cell_h * rows + (40 if title else 0)
    doc = fitz.open()
    page = doc.new_page(width=W, height=H)
    y0 = 0
    if title:
        page.insert_text((30, 26), title, fontsize=13, fontname="hebo")
        y0 = 40
    for i, (path, label) in enumerate(panels):
        r, c = divmod(i, cols)
        x, y = c * cell_w, y0 + r * cell_h
        src = fitz.open(path)
        sp = src[0]
        rect = fitz.Rect(x + 25, y + 30, x + cell_w - 15, y + cell_h - 10)
        page.show_pdf_page(rect, src, 0)
        if label:
            page.insert_text((x + 8, y + 18), label, fontsize=11, fontname="hebo")
    doc.save(out_pdf, garbage=3, deflate=True)
    print("written:", out_pdf)

FK = os.path.join(BASE, "figures_kncc")
FR = os.path.join(BASE, "figures_reframe")

# Figure 1: study design + QC + heterogeneity
compose([
    (f"{FK}/qc_pca_timepoint.pdf", "a"),
    (f"{FK}/qc_pca_plex.pdf", "b"),
    (f"{FK}/techvariance_direct.pdf", "c"),
    (f"{FK}/heterogeneity_v4.pdf", "d"),
], f"{OUT}/Figure_1.pdf", cols=2,
    title="Figure 1. Study design, error-aware QC, and inter-patient heterogeneity")

# Figure 2: two-axis evolution map
compose([
    (f"{FK}/module_coevolution_heatmap.pdf", "a"),
    (f"{FK}/module_evolution_pca.pdf", "b"),
    (f"{FK}/module_jaccard.pdf", "c"),
], f"{OUT}/Figure_2.pdf", cols=2,
    title="Figure 2. Patient-resolved evolution: two dominant, largely orthogonal axes")

# Figure 3: competing-hypothesis adjudication
compose([
    (f"{FK}/global_module_adjudication.pdf", "a"),
    (f"{FK}/composition_sensitivity.pdf", "b"),
    (f"{FK}/m11_sp_vs_prolif_delta.pdf", "c"),
], f"{OUT}/Figure_3.pdf", cols=2,
    title="Figure 3. Competing-hypothesis adjudication of module changes")

# Figure 4: spliceosome worked example + complex controls
compose([
    (f"{FK}/m2_paired_volcano.pdf", "a"),
    (f"{FK}/m3_subcomplex_deltas.pdf", "b"),
    (f"{FK}/complex_controls.pdf", "c"),
    (f"{FK}/glass_second_axis.pdf", "d"),
], f"{OUT}/Figure_4.pdf", cols=2,
    title="Figure 4. Worked example: the spliceosome and matched complex controls")

# Figure 5: phosphoproteome
compose([
    (f"{FK}/m5_raw_vs_corrected.pdf", "a"),
    (f"{FK}/phospho_joint_quadrant.pdf", "b"),
], f"{OUT}/Figure_5.pdf", cols=2,
    title="Figure 5. Phosphoproteome: no abundance-independent changes detected at available power")

# Figure 6: independent validation + spatial + boundaries
compose([
    (f"{FK}/hhs_module_consistency.pdf", "a"),
    (f"{FK}/hhs_pc2_projection.pdf", "b"),
    (f"{FK}/hhs_emt_neuronal.pdf", "c"),
    (f"{FK}/geomx_pc2_by_stain.pdf", "d"),
    (f"{FK}/geomx_pc2_microglia.pdf", "e"),
    (f"{FK}/cptac_matched_validation.pdf", "f"),
], f"{OUT}/Figure_6.pdf", cols=3,
    title="Figure 6. Independent longitudinal validation, spatial context, and boundaries")

print("all figures composed.")
