# Deep cross-file audit: figure refs, legend-file mapping, table refs
import re, os, zipfile

print("=" * 70)
print("DEEP AUDIT: references vs actual files")
print("=" * 70)

ms = open("MANUSCRIPT_GND.md", encoding="utf-8").read()
supp = open("submission_GND/supplementary_materials.md", encoding="utf-8").read()

fails = []

# 1. every "Fig. N" cited in text has a legend AND a file
text_mains = sorted(set(int(x) for x in re.findall(r"Fig\. (\d+)\b", ms)))
print(f"Main figures cited in text: {text_mains}")
for n in text_mains:
    has_legend = f"**Figure {n}." in ms
    has_file = os.path.exists(f"submission_GND/02_Figures/Figure_{n}.jpg")
    if not (has_legend and has_file):
        fails.append(f"Figure {n}: legend={has_legend} file={has_file}")

# 2. every "Fig. SN" cited in text has legend in supplementary AND a pdf
text_supps = sorted(set(int(x) for x in re.findall(r"Fig\. S(\d+)", ms)))
print(f"Supplementary figures cited in text: {text_supps}")
supp_pdfs = os.listdir("submission_GND/Supplementary/figures/pdf_ordered")
for n in text_supps:
    has_legend = f"Figure S{n} |" in supp
    has_file = any(f.startswith(f"S{n:02d}_") for f in supp_pdfs)
    if not (has_legend and has_file):
        fails.append(f"Figure S{n}: legend={has_legend} file={has_file}")

# 3. every legend in supplementary has a pdf
legends = re.findall(r"\*\*Figure (S\d+) \|", supp)
print(f"Legends in supplementary: {len(legends)}")
for n in [int(x[1:]) for x in legends]:
    has_file = any(f.startswith(f"S{n:02d}_") or f.startswith(f"S{n}b" if n == 21 else f"S{n:02d}_") for f in supp_pdfs)
    if n == 21:  # S21b exists too
        has_file = any(f.startswith("S21") for f in supp_pdfs)
    if not has_file:
        fails.append(f"Legend S{n}: no pdf file")

# 4. supplementary tables referenced in text exist in xlsx
names = re.findall(r'name="([^"]+)"', zipfile.ZipFile("submission_GND/04_Supplementary/Supplementary_Tables_GND.xlsx").read("xl/workbook.xml").decode())
text_tables = set(re.findall(r"Supplementary Table (S_\w+|Sthresh|Sperf|Sablation)", ms))
print(f"Tables referenced in text: {text_tables}")
for t in text_tables:
    if t not in names:
        fails.append(f"Table {t}: referenced but no sheet")

# 5. legend-listed tables exist as sheets
legend_tables = re.findall(r"\| (S_\w+) \|", supp)
missing = [t for t in legend_tables if t not in names]
if missing:
    fails.append(f"Legend tables without sheets: {missing}")
extra = [t for t in names if t not in legend_tables]
if extra:
    print(f"NOTE: sheets without legend entries: {extra}")

# 6. main figure count vs journal limit
n_main = len([f for f in os.listdir("submission_GND/02_Figures") if f.endswith(".jpg")])
print(f"Main figures in package: {n_main} (+1 table = {n_main + 1} items; G&D limit 8)")

print("-" * 70)
if fails:
    print(f"FAILED ({len(fails)}):")
    for f in fails:
        print("  -", f)
else:
    print("ALL DEEP CHECKS PASSED")
