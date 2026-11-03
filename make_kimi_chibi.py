# Chibi humanoid Kimi pet sprite — big head, blue hoodie, idle/blink/talk/happy
from PIL import Image, ImageDraw
import os

OUT = r"E:\KIMI\kimi_pet\sprites_chibi"
os.makedirs(OUT, exist_ok=True)

S = 160
SKIN = (255, 225, 200)
HAIR = (60, 70, 90)          # dark navy hair
HOODIE = (43, 109, 255)      # Kimi blue hoodie
HOODIE_D = (30, 85, 210)
DARK = (35, 40, 55)
WHITE = (255, 255, 255)
PINK = (255, 165, 185)
SHOE = (45, 55, 75)

def base():
    img = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    cx = S // 2
    # body (hoodie)
    d.rounded_rectangle([cx-26, 96, cx+26, 148], radius=14, fill=HOODIE)
    d.line([cx, 100, cx, 140], fill=HOODIE_D, width=3)          # zipper
    # hood behind head
    d.pieslice([cx-32, 26, cx+32, 96], start=180, end=360, fill=HOODIE)
    # head
    d.ellipse([cx-30, 30, cx+30, 94], fill=SKIN)
    # hair: top cap + bangs
    d.pieslice([cx-31, 26, cx+31, 80], start=180, end=360, fill=HAIR)
    d.polygon([(cx-31, 52), (cx-31, 66), (cx-20, 58)], fill=HAIR)
    d.polygon([(cx+31, 52), (cx+31, 66), (cx+20, 58)], fill=HAIR)
    d.polygon([(cx-18, 40), (cx-6, 58), (cx-14, 60)], fill=HAIR)
    d.polygon([(cx+18, 40), (cx+6, 58), (cx+14, 60)], fill=HAIR)
    # ahoge (little hair flick, Kimi blue)
    d.line([cx, 28, cx+6, 16], fill=HOODIE, width=3)
    d.line([cx+6, 16, cx+12, 22], fill=HOODIE, width=3)
    # arms
    d.ellipse([cx-36, 100, cx-22, 124], fill=HOODIE)
    d.ellipse([cx+22, 100, cx+36, 124], fill=HOODIE)
    # legs + shoes
    d.rectangle([cx-16, 146, cx-4, 152], fill=DARK)
    d.rectangle([cx+4, 146, cx+16, 152], fill=DARK)
    d.ellipse([cx-19, 148, cx-2, 156], fill=SHOE)
    d.ellipse([cx+2, 148, cx+19, 156], fill=SHOE)
    # cheeks
    d.ellipse([cx-25, 72, cx-16, 79], fill=PINK)
    d.ellipse([cx+16, 72, cx+25, 79], fill=PINK)
    return img, d, cx

def eyes(d, cx, mode="open"):
    if mode == "open":
        for sx in (-14, 14):
            d.ellipse([cx+sx-6, 58, cx+sx+6, 74], fill=WHITE, outline=DARK, width=2)
            d.ellipse([cx+sx-3, 62, cx+sx+3, 72], fill=(70, 130, 255))  # blue iris
            d.ellipse([cx+sx-2, 63, cx+sx+2, 69], fill=DARK)
            d.ellipse([cx+sx-2, 63, cx+sx-0, 65], fill=WHITE)
    elif mode == "blink":
        for sx in (-14, 14):
            d.line([cx+sx-5, 66, cx+sx+5, 66], fill=DARK, width=3)
    elif mode == "happy":
        for sx in (-14, 14):
            d.arc([cx+sx-6, 58, cx+sx+6, 74], start=180, end=360, fill=DARK, width=3)

def mouth(d, cx, kind):
    if kind == "smile":
        d.arc([cx-5, 72, cx+5, 82], start=0, end=180, fill=DARK, width=2)
    elif kind == "talk":
        d.ellipse([cx-5, 74, cx+5, 84], fill=DARK)
        d.ellipse([cx-3, 77, cx+3, 82], fill=PINK)
    elif kind == "happy":
        d.pieslice([cx-6, 70, cx+6, 84], start=0, end=180, fill=DARK)

img, d, cx = base(); eyes(d, cx, "open"); mouth(d, cx, "smile");  img.save(f"{OUT}/idle.png")
img, d, cx = base(); eyes(d, cx, "blink"); mouth(d, cx, "smile"); img.save(f"{OUT}/blink.png")
img, d, cx = base(); eyes(d, cx, "open"); mouth(d, cx, "talk");   img.save(f"{OUT}/talk.png")
img, d, cx = base(); eyes(d, cx, "happy"); mouth(d, cx, "happy"); img.save(f"{OUT}/happy.png")

# preview strip
imgs = [Image.open(f"{OUT}/{n}.png") for n in ["idle", "blink", "talk", "happy"]]
w = sum(i.width for i in imgs); h = max(i.height for i in imgs)
strip = Image.new("RGBA", (w, h), (255, 255, 255, 255))
x = 0
for i in imgs: strip.paste(i, (x, 0), i); x += i.width
strip.save(f"{OUT}/_preview.png")
print("chibi sprites done")
