# Final cross-file consistency audit (manuscript vs results vs supplementary)
import re, os

print("=" * 70)
print("FINAL CROSS-FILE AUDIT")
print("=" * 70)

ms = open("MANUSCRIPT_GND.md", encoding="utf-8").read()
val = open("results/sixgene_validation_summary.csv").read()
perf = open("results/signature_performance_full.txt").read()

checks = []
def chk(name, cond, detail=""):
    checks.append((name, cond, detail))

# 1. six-gene numbers: manuscript KM p must match validation CSV km_p
chk("CGGA-693 KM p=3.4e-19 in text & CSV",
    "3.4 × 10⁻¹⁹" in ms and "3.40239202388763e-19" in val)
chk("CGGA-325 KM p=7.7e-18 in text & CSV",
    "7.7 × 10⁻¹⁸" in ms and "7.73411337353704e-18" in val)
chk("GLASS HR=2.14 p=2e-06 in text & glass_validation",
    "HR = 2.14" in ms and "2.138498" in open("results/glass_validation.txt").read())

# 2. stale value must NOT appear anywhere
chk("No stale Cox p=0.948 for CGGA-325", "0.948281426547432" not in val)

# 3. Virtual KO eliminated
chk("No 'Virtual KO' in manuscript", "Virtual KO" not in ms and "Virtual knockout" not in ms)

# 4. meta-analysis wording
chk("Random-effects meta (not fixed)", "random-effects meta-analysis" in ms)
chk("I2 reported", "91.4%" in ms)

# 5. maxRank wording
chk("maxRank=1500 = 7.5% universe", "7.5% of the ~20,000-gene universe" in ms)

# 6. bootstrap 1000
chk("Bootstrap 1,000 in text", "1,000 iterations" in ms)

# 7. GLASS units
chk("GLASS 213 pairs/173 patients in text", "213 pairs from 173 independent patients" in ms or "213 pairs/173" in ms)
chk("179 pairs/148 patients in text", "179 pairs from 148 patients" in ms)

# 8. 107/111 vs 110/111
chk("107 sig + 110 positive logFC distinguished", "110 of the 111 proteins showed positive log" in ms)

# 9. supplementary sheets
import zipfile
names = re.findall(r'name="([^"]+)"', zipfile.ZipFile("submission_GND/04_Supplementary/Supplementary_Tables_GND.xlsx").read("xl/workbook.xml").decode())
chk("S_ablation sheet exists", "S_ablation" in names)
chk("S_Table2_targets exists", "S_Table2_targets" in names)
chk("S_glass_design exists", "S_glass_design" in names)
chk("S_splicing_sig4387 exists", "S_splicing_sig4387" in names)

# 10. supp front title
supp = open("submission_GND/supplementary_materials.md", encoding="utf-8").read()
chk("Supp front-page new title", "marks a writer-dominant" in supp)
chk("S24-S33 legends present", "Figure S33" in supp and "Figure S24" in supp)

ok = sum(1 for _, c, _ in checks if c)
for name, cond, detail in checks:
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
print("-" * 70)
print(f"{ok}/{len(checks)} checks passed")
