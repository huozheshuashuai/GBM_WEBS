# =============================================================================
# Kimi desktop pet sprite generator
# A round Kimi-blue blob with big friendly eyes — the Kimi you know
# Frames: idle, blink, talk, happy
# =============================================================================
from PIL import Image, ImageDraw
import os

OUT = r"E:\KIMI\kimi_pet\sprites"
os.makedirs(OUT, exist_ok=True)

S = 128
KIMI_BLUE = (43, 109, 255)      # Kimi brand blue
KIMI_LIGHT = (120, 170, 255)
DARK = (30, 40, 60)
WHITE = (255, 255, 255)
PINK = (255, 170, 190)

def base_frame():
    img = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    # body: round blob, Kimi blue with lighter belly shading
    d.ellipse([20, 26, 108, 118], fill=KIMI_BLUE)
    d.ellipse([30, 36, 98, 108], fill=KIMI_LIGHT)     # inner lighter area
    d.ellipse([40, 74, 88, 106], fill=(230, 240, 255))  # belly
    # little feet
    d.ellipse([40, 112, 58, 122], fill=DARK)
    d.ellipse([70, 112, 88, 122], fill=DARK)
    # small tuft antenna on head (Kimi's little sprout)
    d.line([64, 26, 64, 14], fill=DARK, width=3)
    d.ellipse([58, 6, 70, 16], fill=KIMI_BLUE, outline=DARK)
    # cheeks
    d.ellipse([32, 72, 44, 82], fill=PINK)
    d.ellipse([84, 72, 96, 82], fill=PINK)
    return img, d

def eyes_open(d, happy=False):
    if happy:
        d.arc([40, 48, 56, 66], start=180, end=360, fill=DARK, width=3)
        d.arc([72, 48, 88, 66], start=180, end=360, fill=DARK, width=3)
    else:
        d.ellipse([42, 50, 56, 68], fill=WHITE, outline=DARK, width=2)
        d.ellipse([72, 50, 86, 68], fill=WHITE, outline=DARK, width=2)
        d.ellipse([46, 54, 54, 66], fill=DARK)
        d.ellipse([76, 54, 84, 66], fill=DARK)
        d.ellipse([48, 55, 51, 59], fill=WHITE)
        d.ellipse([78, 55, 81, 59], fill=WHITE)

def mouth(d, kind):
    if kind == "smile":
        d.arc([56, 66, 72, 78], start=0, end=180, fill=DARK, width=2)
    elif kind == "talk":
        d.ellipse([57, 68, 71, 80], fill=DARK)
        d.ellipse([59, 71, 69, 78], fill=PINK)
    elif kind == "happy":
        d.pieslice([55, 64, 73, 80], start=0, end=180, fill=DARK)

img, d = base_frame(); eyes_open(d); mouth(d, "smile");  img.save(os.path.join(OUT, "idle.png"))
img, d = base_frame()
d.line([44, 60, 54, 60], fill=DARK, width=3); d.line([74, 60, 84, 60], fill=DARK, width=3)
mouth(d, "smile"); img.save(os.path.join(OUT, "blink.png"))
img, d = base_frame(); eyes_open(d); mouth(d, "talk");   img.save(os.path.join(OUT, "talk.png"))
img, d = base_frame(); eyes_open(d, happy=True); mouth(d, "happy"); img.save(os.path.join(OUT, "happy.png"))

print("kimi sprites done:", os.listdir(OUT))
