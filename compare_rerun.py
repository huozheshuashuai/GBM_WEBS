# Compare snapshot vs re-run results: extract key numbers and check consistency
import re, os

KEY_FILES = [
    "cptac_summary.txt", "acetyl_summary.txt", "phospho_summary.txt",
    "splicing_summary.txt", "splicing_lgg_summary.txt", "splicing_tme_summary.txt",
    "virtual_ko_abs.txt", "ptm_abundance_correction.txt", "uncoupling_rigor.txt",
    "incremental_value.txt", "glass_validation.txt", "glass_recurrence_time.txt",
    "glass_firstpair_sensitivity.txt", "gse83300_validation.txt",
    "threshold_sensitivity.txt", "signature_performance_full.txt",
    "bootstrap_stability_summary.txt",
]

def numbers(path):
    if not os.path.exists(path):
        return None
    txt = open(path, encoding="utf-8", errors="replace").read()
    return re.findall(r"-?\d+\.?\d*(?:e[+-]?\d+)?", txt, flags=re.I)

print(f"{'file':45s} {'match':>6s} {'n_old':>6s} {'n_new':>6s}")
print("-" * 70)
for f in KEY_FILES:
    old = numbers(f"results_snapshot_round1/{f}")
    new = numbers(f"results/{f}")
    if old is None or new is None:
        print(f"{f:45s} {'MISSING':>6s}")
        continue
    n_same = sum(1 for a, b in zip(old, new) if a == b)
    pct = 100 * n_same / max(len(old), 1)
    flag = "OK" if pct > 95 else ("~" if pct > 80 else "DIFF")
    print(f"{f:45s} {pct:5.0f}% {len(old):>6d} {len(new):>6d}  {flag}")
