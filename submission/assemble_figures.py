# Final assembly: composite Figure_1..Figure_12 in jpg/tif/pdf
import fitz, os, shutil
from PIL import Image

BASE = r"E:\KIMI\GBM_WEBS_Project"
SUB = os.path.join(BASE, "submission")
FIGDIR = os.path.join(SUB, "Figures")
shutil.rmtree(FIGDIR, ignore_errors=True)
for fmt in ("jpg", "tif", "pdf"):
    os.makedirs(os.path.join(FIGDIR, fmt), exist_ok=True)

def merge_pdfs(paths, out_pdf):
    doc = fitz.open()
    for p in paths:
        d = fitz.open(p); doc.insert_pdf(d); d.close()
    doc.save(out_pdf); doc.close()

def rasterize(pdf_path, outdir, stem, dpi=300):
    doc = fitz.open(pdf_path)
    for i, page in enumerate(doc):
        sfx = stem if len(doc) == 1 else f"{stem}_p{i+1}"
        pix = page.get_pixmap(dpi=dpi)
        img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
        img.save(os.path.join(outdir, "jpg", sfx + ".jpg"), quality=95)
        img.save(os.path.join(outdir, "tif", sfx + ".tif"), compression="tiff_lzw")
    doc.close()

def rasterize_png(png_path, outdir, stem):
    img = Image.open(png_path).convert("RGB")
    img.save(os.path.join(outdir, "jpg", stem + ".jpg"), quality=95)
    img.save(os.path.join(outdir, "tif", stem + ".tif"), compression="tiff_lzw")
    shutil.copy(png_path, os.path.join(outdir, "pdf", stem + ".png"))

F = os.path.join(BASE, "figures")
FV2 = os.path.join(BASE, "figures_v2")
FP = os.path.join(BASE, "figures_paper")
FF = os.path.join(BASE, "figures_fixed")
ORIG = r"E:\KIMI\extracted"

figure_map = {
    "Figure_1":  ("png",  os.path.join(ORIG, "Figure1_StudyDesign_v16_modified.png")),
    "Figure_2":  ("pdf",  [os.path.join(FP, "Fig2_WEBS_landscape.pdf")]),
    "Figure_3":  ("pdf",  [os.path.join(FP, "Fig3_Neftel_WEBS.pdf")]),
    "Figure_4":  ("pdf",  [os.path.join(FP, "Fig4_hypoxia_decoupling.pdf")]),
    "Figure_5":  ("pdf",  [os.path.join(F, "fig_pseudotime.pdf")]),
    "Figure_6":  ("pdf",  [os.path.join(FF, "fig_cellchat_compare.pdf"),
                           os.path.join(FF, "fig_cellchat_ranknet.pdf"),
                           os.path.join(FF, "fig_cellchat_diffnet.pdf")]),
    "Figure_7":  ("pdf",  [os.path.join(FF, f"fig_spatial_{s}.pdf") for s in
                           ["GBM1","GBM2","GBM3","GBM4","GBM5_1","GBM5_2"]]),
    "Figure_8":  ("pdf",  [os.path.join(FV2, "Fig8_snATAC_promoter.pdf")]),
    "Figure_9":  ("pdf",  [os.path.join(FP, "Fig9_functional.pdf")]),
    "Figure_10": ("pdf",  [os.path.join(FV2, "Fig10_external_validation.pdf"),
                           os.path.join(FF, "fig_signature_TCGA_GBM.pdf"),
                           os.path.join(FF, "fig_signature_cgga693.pdf"),
                           os.path.join(FF, "fig_signature_cgga325.pdf"),
                           os.path.join(F, "fig_signature_lasso.pdf")]),
    "Figure_11": ("pdf",  [os.path.join(F, "fig_depmap_scatter.pdf"),
                           os.path.join(F, "fig_gdsc2_scatter.pdf"),
                           os.path.join(F, "fig_cmap_drugs.pdf")]),
    "Figure_12": ("png",  os.path.join(ORIG, "Figure12_WorkingModel_v16.png")),
}

tmp = os.path.join(SUB, "_tmp"); os.makedirs(tmp, exist_ok=True)
for name, (kind, src) in figure_map.items():
    if kind == "png":
        rasterize_png(src, FIGDIR, name)
    else:
        out_pdf = os.path.join(tmp, name + ".pdf")
        merge_pdfs(src, out_pdf)
        shutil.copy(out_pdf, os.path.join(FIGDIR, "pdf", name + ".pdf"))
        rasterize(out_pdf, FIGDIR, name)
    print(name, "done")

shutil.rmtree(tmp)
print("ALL FIGURES ASSEMBLED")
