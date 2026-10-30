# Composite multi-part figures into single-canvas images (jpg/tif/pdf)
import fitz, os
from PIL import Image

BASE = r"E:\KIMI\GBM_WEBS_Project"
SUB = os.path.join(BASE, "submission")
FIGDIR = os.path.join(SUB, "Figures")

def render_pdf_pages(pdf_path, dpi=200):
    doc = fitz.open(pdf_path)
    imgs = []
    for page in doc:
        pix = page.get_pixmap(dpi=dpi)
        imgs.append(Image.frombytes("RGB", (pix.width, pix.height), pix.samples))
    doc.close()
    return imgs

def make_grid(images, ncol, out_stem, pad=30, label=True):
    # uniform cell size
    cw = max(im.width for im in images)
    ch = max(im.height for im in images)
    nrow = (len(images) + ncol - 1) // ncol
    W = ncol * cw + (ncol + 1) * pad
    H = nrow * ch + (nrow + 1) * pad
    canvas = Image.new("RGB", (W, H), "white")
    labels = "abcdefgh"
    from PIL import ImageDraw
    draw = ImageDraw.Draw(canvas)
    for i, im in enumerate(images):
        r, c = divmod(i, ncol)
        x = pad + c * (cw + pad)
        y = pad + r * (ch + pad)
        canvas.paste(im, (x, y))
    canvas.save(out_stem + ".jpg", quality=95)
    canvas.save(out_stem + ".tif", compression="tiff_lzw")
    # single-page PDF from image
    img_pdf = fitz.open()
    rect = fitz.Rect(0, 0, canvas.width * 72 / 200, canvas.height * 72 / 200)
    page = img_pdf.new_page(width=rect.width, height=rect.height)
    import io
    buf = io.BytesIO()
    canvas.save(buf, format="PNG")
    page.insert_image(rect, stream=buf.getvalue())
    img_pdf.save(out_stem + ".pdf")
    img_pdf.close()

composites = {
    "Figure_6":  (["fig_cellchat_compare.pdf", "fig_cellchat_ranknet.pdf",
                   "fig_cellchat_diffnet.pdf"], "figures_fixed", 3),
    "Figure_7":  ([f"fig_spatial_{s}.pdf" for s in
                   ["GBM1","GBM2","GBM3","GBM4","GBM5_1","GBM5_2"]], "figures_fixed", 3),
    "Figure_10": (["Fig10_external_validation.pdf", "fig_signature_TCGA_GBM.pdf",
                   "fig_signature_cgga693.pdf", "fig_signature_cgga325.pdf",
                   "fig_signature_lasso.pdf"], None, 3),  # mixed dirs
    "Figure_11": (["fig_depmap_scatter.pdf", "fig_gdsc2_scatter.pdf",
                   "fig_cmap_drugs.pdf"], "figures", 3),
}

srcdir = {"figures_fixed": os.path.join(BASE, "figures_fixed"),
          "figures": os.path.join(BASE, "figures")}
fv2 = os.path.join(BASE, "figures_v2")

for name, (files, default_dir, ncol) in composites.items():
    images = []
    for f in files:
        if default_dir and os.path.exists(os.path.join(srcdir[default_dir], f)):
            p = os.path.join(srcdir[default_dir], f)
        elif os.path.exists(os.path.join(fv2, f)):
            p = os.path.join(fv2, f)
        elif os.path.exists(os.path.join(srcdir["figures_fixed"], f)):
            p = os.path.join(srcdir["figures_fixed"], f)
        else:
            p = os.path.join(BASE, "figures", f)
        images.extend(render_pdf_pages(p))
    make_grid(images, ncol, os.path.join(FIGDIR, "jpg", name))
    # move pdf into pdf/ and tif into tif/
    os.replace(os.path.join(FIGDIR, "jpg", name + ".pdf"), os.path.join(FIGDIR, "pdf", name + ".pdf"))
    os.replace(os.path.join(FIGDIR, "jpg", name + ".tif"), os.path.join(FIGDIR, "tif", name + ".tif"))
    print(name, "composited:", len(images), "panels")

# remove stale multi-page splits for composited figures
for name in composites:
    for fmt in ("jpg", "tif"):
        for f in os.listdir(os.path.join(FIGDIR, fmt)):
            if f.startswith(name + "_p"):
                os.remove(os.path.join(FIGDIR, fmt, f))
print("done")
