# Build merged Supplementary PDF: text (methods+legends) + figure pages
import fitz, os, re

BASE = r"E:\KIMI\GBM_WEBS_Project"
SUB = os.path.join(BASE, "submission")
MD = os.path.join(SUB, "supplementary_materials.md")
FIGDIR = os.path.join(SUB, "Supplementary", "figures", "pdf")
OUT = os.path.join(SUB, "Supplementary", "Supplementary_Materials_merged.pdf")

def md_to_lines(path):
    """Convert simple markdown to (text, style) lines."""
    lines = []
    for raw in open(path, encoding="utf-8"):
        line = raw.rstrip("\n")
        if line.startswith("# "):
            lines.append((line[2:], "h1"))
        elif line.startswith("## "):
            lines.append((line[3:], "h2"))
        elif line.startswith("**") and line.endswith("**"):
            lines.append((line.strip("*"), "b"))
        elif line.strip() == "---":
            lines.append(("", "sep"))
        else:
            # strip md bold/italic markers and table pipes to plain text
            line = re.sub(r"\*\*(.+?)\*\*", r"\1", line)
            lines.append((line, "p"))
    return lines

def add_text_pages(doc, lines):
    W, H = fitz.paper_size("a4")
    margin = 56
    page = doc.new_page(width=W, height=H)
    y = margin
    for text, style in lines:
        if style == "sep":
            y += 10
            continue
        size = {"h1": 16, "h2": 13, "b": 11, "p": 10}[style]
        font = "hebo" if style in ("h1", "h2", "b") else "helv"
        if text == "":
            y += size * 0.8
            continue
        rect = fitz.Rect(margin, y, W - margin, y + 500)
        rc = page.insert_textbox(rect, text, fontsize=size, fontname=font)
        est_h = max(size * 1.5, -rc)
        if rc < 0 and False:
            pass
        # estimate wrapped height
        import math
        chars_per_line = int((W - 2 * margin) / (size * 0.5))
        n_lines = max(1, math.ceil(len(text) / chars_per_line))
        est_h = n_lines * size * 1.35 + 4
        if y + est_h > H - margin:
            page = doc.new_page(width=W, height=H)
            y = margin
            rect = fitz.Rect(margin, y, W - margin, y + 500)
            page.insert_textbox(rect, text, fontsize=size, fontname=font)
        y += est_h
    return doc

doc = fitz.open()
lines = md_to_lines(MD)
doc = add_text_pages(doc, lines)

# append supplementary figure pages
figs = sorted([f for f in os.listdir(FIGDIR) if f.endswith(".pdf")])
for f in figs:
    src = fitz.open(os.path.join(FIGDIR, f))
    doc.insert_pdf(src)
    src.close()

doc.save(OUT)
print(f"Merged PDF: {len(doc)} pages ({len(figs)} figure files) -> {OUT}")
