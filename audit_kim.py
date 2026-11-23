import re, html, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

t = open(r"E:/KIMI/kim2024.html", encoding="utf-8", errors="replace").read()
t = html.unescape(t)
t = re.sub(r"<script.*?</script>", " ", t, flags=re.S)
t = re.sub(r"<style.*?</style>", " ", t, flags=re.S)
text = re.sub(r"<[^>]+>", " ", t)
text = re.sub(r"\s+", " ", text)

# Results section boundaries
i0 = text.find("The molecular trajectory of longitudinal glioblastoma")
i1 = text.find("Discussion", i0)
res = text[i0:i1]
print("results length:", len(res))

# scan for key methods/analyses
keys = ["phospho", "KSEA", "kinase", "trajectory", "heterogen", "predict", "WGCNA", "MOFA", "multi-omics factor",
        "cluster", "subtype", "splice", "splicing", "metaboli", "OXPHOS", "immune", "microenvironment",
        "hypermutat", "MMR", "MGMT", "temozolomide", "TMZ", "survival", "prognos", "PFS", "overall survival",
        "patient-specific", "per-patient", "intertumor", "inter-tumor", "variance", "acetylat", "ubiquitin"]
print("\n=== keyword hits in Results (with context) ===")
for k in keys:
    hits = [m.start() for m in re.finditer(k, res, re.I)]
    print(f"\n### {k}: {len(hits)} hits")
    for h in hits[:3]:
        print("   ...", res[max(0,h-140):h+160].strip()[:320])
