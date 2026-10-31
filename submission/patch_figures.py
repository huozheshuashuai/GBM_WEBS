# Patch Figure 1 & 12: update numbers and Chinese->English labels
import os
from PIL import Image, ImageDraw, ImageFont

FONT = r"C:\Windows\Fonts\arial.ttf"
FONT_B = r"C:\Windows\Fonts\arialbd.ttf"
ORIG = r"E:\KIMI\extracted"
OUT = r"E:\KIMI\GBM_WEBS_Project\figures_final"
os.makedirs(OUT, exist_ok=True)

def bg_sample(img, box, margin=8):
    """sample background color from a strip above the box"""
    x0, y0, x1, y1 = box
    strip = img.crop((x0, max(0, y0 - margin - 4), x1, max(0, y0 - 4)))
    px = list(strip.getdata())
    px.sort()
    return px[len(px) // 2]

def patch(img, box, new_text, font_path, font_size, text_color, align="center",
          bg_color=None, line_spacing=4):
    if bg_color is None:
        bg_color = bg_sample(img, box)
    draw = ImageDraw.Draw(img)
    draw.rectangle(box, fill=bg_color)
    font = ImageFont.truetype(font_path, font_size)
    lines = new_text.split("\n")
    heights = [draw.textbbox((0, 0), l, font=font)[3] for l in lines]
    total_h = sum(heights) + line_spacing * (len(lines) - 1)
    y = box[1] + (box[3] - box[1] - total_h) / 2
    for i, l in enumerate(lines):
        w = draw.textbbox((0, 0), l, font=font)[2]
        if align == "center":
            x = box[0] + (box[2] - box[0] - w) / 2
        else:
            x = box[0] + 6
        draw.text((x, y), l, font=font, fill=text_color)
        y += heights[i] + line_spacing
    return img

# =============== Figure 12 ===============
img12 = Image.open(os.path.join(ORIG, "Figure12_WorkingModel_v16.png")).convert("RGB")

# 1. LLS Chinese label -> English
patch(img12, (150, 115, 700, 190), "LLS  low-balance state", FONT_B, 34, (26, 32, 44))
# 2. LHS Chinese label -> English
patch(img12, (1380, 115, (1980), 190), "LHS  high-balance state", FONT_B, 34, (26, 32, 44))
# 3. spliceosome numbers
patch(img12, (1560, 285, 2040, 340), "107/111 spliceosome proteins\nup-regulated (p = 3.4e-37)",
      FONT_B, 26, (180, 40, 40))
# 4. vorinostat rho -> q
patch(img12, (245, 1048, 820, 1115), "Vorinostat/SAHA\nHigh-WEBS sensitivity q = 1.1e-4",
      FONT, 24, (20, 20, 20), align="left")
# 5. MLLT1/ENL -> remove (cover with box purple)
patch(img12, (1908, 430, 2044, 472), "", FONT, 20, (255, 255, 255),
      bg_color=(147, 112, 183))

img12.save(os.path.join(OUT, "Figure12_WorkingModel_final.png"))

# =============== Figure 1 ===============
img1 = Image.open(os.path.join(ORIG, "Figure1_StudyDesign_v16_modified.png")).convert("RGB")

# 1. cell count
patch(img1, (95, 52, 330, 80), "17 patients · 35,823 cells", FONT, 22, (60, 50, 80),
      align="left")
# 2. CPTAC numbers (two lines inside the orange box)
patch(img1, (1445, 822, 1850, 862), "107 / 111 proteins up-regulated", FONT_B, 21,
      (60, 40, 110), align="left")
patch(img1, (1445, 860, 1880, 900), "p = 3.4e-37  (top enriched pathway)", FONT_B, 21,
      (60, 40, 110), align="left")

img1.save(os.path.join(OUT, "Figure1_StudyDesign_final.png"))

print("Both figures patched ->", OUT)
