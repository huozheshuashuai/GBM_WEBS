# =============================================================================
# 小WE (Wee) — desktop pet sprite generator
# A round purple-teal creature with a tiny balance scale on its head
# Frames: idle, blink, talk, happy
# =============================================================================
from PIL import Image, ImageDraw
import os

OUT = r"E:\KIMI\deskpet_wee\sprites"
os.makedirs(OUT, exist_ok=True)

S = 128  # canvas size
PURPLE = (155, 124, 200)      # writers side
TEAL = (95, 180, 175)         # erasers side
DARK = (40, 40, 50)
WHITE = (255, 255, 255)
PINK = (255, 170, 180)

def base_frame():
    img = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    # body: two-tone blob
    d.ellipse([18, 34, 110, 116], fill=PURPLE)          # full body
    d.pieslice([18, 34, 110, 116], start=-90, end=90, fill=TEAL)  # right half teal
    # belly
    d.ellipse([42, 78, 86, 110], fill=(245, 240, 250))
    # feet
    d.ellipse([38, 110, 58, 122], fill=DARK)
    d.ellipse([70, 110, 90, 122], fill=DARK)
    # balance scale on head: stick + beam + two pans
    d.line([64, 10, 64, 36], fill=DARK, width=3)        # stick
    d.line([42, 14, 86, 14], fill=DARK, width=3)        # beam
    d.line([42, 14, 38, 26], fill=DARK, width=2)        # left string
    d.line([86, 14, 90, 26], fill=DARK, width=2)        # right string
    d.pieslice([32, 22, 46, 36], start=0, end=180, fill=PURPLE, outline=DARK)   # left pan
    d.pieslice([82, 22, 96, 36], start=0, end=180, fill=TEAL, outline=DARK)     # right pan
    d.ellipse([61, 8, 67, 14], fill=GOLD if False else (255, 210, 90), outline=DARK)  # top knob
    # cheeks
    d.ellipse([28, 74, 42, 84], fill=PINK)
    d.ellipse([86, 74, 100, 84], fill=PINK)
    return img, d

def eyes_open(d, happy=False):
    if happy:
        # happy arc eyes ^_^
        d.arc([38, 48, 56, 68], start=180, end=360, fill=DARK, width=3)
        d.arc([72, 48, 90, 68], start=180, end=360, fill=DARK, width=3)
    else:
        d.ellipse([40, 50, 56, 70], fill=WHITE, outline=DARK, width=2)
        d.ellipse([72, 50, 88, 70], fill=WHITE, outline=DARK, width=2)
        d.ellipse([45, 55, 53, 66], fill=DARK)
        d.ellipse([77, 55, 85, 66], fill=DARK)
        d.ellipse([47, 56, 50, 60], fill=WHITE)  # sparkle
        d.ellipse([79, 56, 82, 60], fill=WHITE)

def mouth(d, kind):
    if kind == "smile":
        d.arc([56, 66, 72, 80], start=0, end=180, fill=DARK, width=2)
    elif kind == "talk":
        d.ellipse([56, 68, 72, 82], fill=DARK)
        d.ellipse([58, 72, 70, 80], fill=PINK)
    elif kind == "happy":
        d.pieslice([54, 64, 74, 82], start=0, end=180, fill=DARK)

# idle
img, d = base_frame(); eyes_open(d); mouth(d, "smile")
img.save(os.path.join(OUT, "idle.png"))
# blink
img, d = base_frame()
d.line([42, 60, 54, 60], fill=DARK, width=3)
d.line([74, 60, 86, 60], fill=DARK, width=3)
mouth(d, "smile")
img.save(os.path.join(OUT, "blink.png"))
# talk
img, d = base_frame(); eyes_open(d); mouth(d, "talk")
img.save(os.path.join(OUT, "talk.png"))
# happy
img, d = base_frame(); eyes_open(d, happy=True); mouth(d, "happy")
img.save(os.path.join(OUT, "happy.png"))

print("sprites done:", os.listdir(OUT))
