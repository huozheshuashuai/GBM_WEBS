import csv, sys, io
from collections import Counter
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

rows = list(csv.DictReader(open(r"E:/KIMI/prjna_runinfo.csv", encoding="utf-8")))
print("total runs:", len(rows))
print(Counter(r["LibraryStrategy"] for r in rows))

rna = [r for r in rows if r["LibraryStrategy"] == "RNA-Seq"]
wxs = [r for r in rows if r["LibraryStrategy"] == "WXS"]
print("\nRNA-Seq runs:", len(rna), "| total GB:", round(sum(float(r["size_MB"]) for r in rna) / 1024, 1))
print("WXS runs:", len(wxs), "| total GB:", round(sum(float(r["size_MB"]) for r in wxs) / 1024, 1))

# sample names like GBM335 / CNU_P10B; check RNA sample naming + layout
print("\nRNA sample name examples:", [r["SampleName"] for r in rna[:10]])
print("avg MB per RNA run:", round(sum(float(r["size_MB"]) for r in rna) / len(rna)))
print("layout:", Counter(r["LibraryLayout"] for r in rna))

# check T1/T2 info availability in SampleName
import re
tp = Counter()
for r in rna:
    m = re.search(r"(T[12]|_P|_R|Primary|Recur)", r["SampleName"] or "")
    tp[r["SampleName"][-2:] if r["SampleName"] else "?"] += 1
print("sample-name suffix counts:", dict(list(tp.items())[:15]))

# smallest RNA runs first for test download
rna_sorted = sorted(rna, key=lambda r: float(r["size_MB"]))
print("\nsmallest RNA run:", rna_sorted[0]["Run"], rna_sorted[0]["size_MB"], "MB", rna_sorted[0]["SampleName"])
print("has download_path:", rna_sorted[0]["download_path"][:80])
