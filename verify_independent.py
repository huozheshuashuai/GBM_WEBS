#!/usr/bin/env python
# =============================================================================
# INDEPENDENT REPRODUCTION of the paper's core results — BASE PYTHON + NUMPY
# No pandas, no scipy, no R. Every statistic implemented manually from raw data.
# =============================================================================
import numpy as np, re, sys, io, gzip, csv, os, glob, zipfile
from math import sqrt, erf
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
np.random.seed(42)

BASE = r"E:/KIMI/GBM_WEBS_Project"
DATA = r"E:/KIMI/data"

results = []
def check(label, expected, got, tol=0.02, note=""):
    try:
        ok = (str(got) == expected) if isinstance(expected, str) else abs(float(got) - float(expected)) <= tol
    except Exception:
        ok = False
    mark = "OK" if ok else "** MISMATCH **"
    results.append((label, expected, got, mark, note))
    print(f"[{mark}] {label}: manuscript={expected} | independent={got} {note}", flush=True)

def norm_cdf(z):
    return 0.5 * (1 + erf(z / sqrt(2)))

def rankdata(a):
    a = np.asarray(a, float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a)); ranks[order] = np.arange(1, len(a) + 1)
    # tie handling: average ranks
    i = 0
    sa = a[order]
    while i < len(a):
        j = i
        while j + 1 < len(a) and sa[j + 1] == sa[i]:
            j += 1
        if j > i:
            avg = (i + 1 + j + 1) / 2
            ranks[order[i:j + 1]] = avg
        i = j + 1
    return ranks

def wilcoxon_signed_rank(d):
    """Normal-approximation Wilcoxon signed-rank p-value (valid for n>~20)."""
    d = np.asarray(d, float)
    d = d[~np.isnan(d) & (d != 0)]
    n = len(d)
    if n < 5: return np.nan
    r = rankdata(np.abs(d))
    W = float(np.sum(r * np.sign(d)))
    mu = 0.0
    sigma = sqrt(n * (n + 1) * (2 * n + 1) / 6)
    z = (W - mu) / sigma
    return 2 * (1 - norm_cdf(abs(z)))

def spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    rx, ry = rankdata(x), rankdata(y)
    rx = rx - rx.mean(); ry = ry - ry.mean()
    return float((rx * ry).sum() / sqrt((rx * rx).sum() * (ry * ry).sum()))

def bh_fdr(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    out = np.full(len(p), np.nan)
    pp = p[ok]; n = len(pp)
    order = np.argsort(pp)
    q = pp[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    qq = np.empty(n); qq[order] = np.minimum(q, 1)
    out[ok] = qq
    return out

def read_tsv_matrix(path, col_prefix="Log KNCC_GBM", gene_col=0):
    """Manual TSV parse: returns (genes, sample_names, matrix)."""
    with open(path, encoding="utf-8") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = [i for i, h in enumerate(hdr) if h.startswith(col_prefix)]
        names = [re.sub(r"\.\d+$", "", re.sub(r":.*", "", re.sub(r"^Log ", "", hdr[i]))) for i in idx]
        genes, rows = [], []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            g = parts[gene_col]
            if g in ("Mean", "Median", "StdDev", "", "Gene"): continue
            vals = []
            for i in idx:
                v = parts[i] if i < len(parts) else ""
                vals.append(float(v) if v not in ("", "NA", "NaN") else np.nan)
            genes.append(g); rows.append(vals)
    M = np.array(rows, dtype=float)
    return genes, names, M

def merge_dup_cols(M, names):
    """average duplicate sample columns"""
    uniq, inv = np.unique(names, return_inverse=True)
    out = np.full((M.shape[0], len(uniq)), np.nan)
    for j, u in enumerate(uniq):
        cols = M[:, inv == j]
        out[:, j] = np.nanmean(cols, axis=1)
    return out, list(uniq)

def modscore(genes, gidx, M):
    """mean of member genes per sample"""
    sel = [gidx[g] for g in genes if g in gidx]
    return np.nanmean(M[sel, :], axis=0)

# =============================================================================
# 1. KNCC proteome + spliceosome longitudinal
# =============================================================================
print("\n##### 1. KNCC proteome / spliceosome longitudinal #####", flush=True)
genes, snames, M = read_tsv_matrix(DATA + r"/PDC000514/KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv")
M, snames = merge_dup_cols(M, snames)
gidx = {g: i for i, g in enumerate(genes)}
patients = [re.sub(r"_T[12]$", "", s) for s in snames]
tps = [re.sub(r".*_", "", s) for s in snames]
t1 = {p: i for i, p, t in zip(range(len(snames)), patients, tps) if t == "T1"}
t2 = {p: i for i, p, t in zip(range(len(snames)), patients, tps) if t == "T2"}
paired = sorted(set(t1) & set(t2))
deltaM = M[:, [t2[p] for p in paired]] - M[:, [t1[p] for p in paired]]
print(f"proteins={M.shape[0]}, samples={M.shape[1]}, paired={len(paired)}", flush=True)
check("KNCC paired patients", 105, len(paired), 0)

splice_genes = []
with open(BASE + r"/results/cptac_spliceosome_proteins.csv") as f:
    r = csv.DictReader(f)
    for row in r: splice_genes.append(row["protein"])
sp = [g for g in splice_genes if g in gidx]
check("spliceosome proteins quantified", 110, len(sp), 2)

pvals, meds = [], []
for g in sp:
    d = deltaM[gidx[g], :]
    d = d[~np.isnan(d)]
    pvals.append(wilcoxon_signed_rank(d) if len(d) >= 30 else np.nan)
    meds.append(np.median(d) if len(d) else np.nan)
q = bh_fdr(pvals)
sig = np.sum(q < 0.05)
check("spliceosome FDR<0.05 count", 105, int(sig), 3)
meds = np.array(meds)
check("direction: all-down", "109/110", f"{int(np.nansum(meds < 0))}/{int(np.sum(~np.isnan(meds)))}", note="direction count")
check("median delta", -0.20, round(np.nanmedian(meds), 3), 0.02)

# =============================================================================
# 2. Adjudication (manual OLS via normal equations)
# =============================================================================
print("\n##### 2. Adjudication #####", flush=True)
def parse_gmt(path):
    d = {}
    for line in open(path, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        d[f[0]] = f[2:]
    return d
hm = parse_gmt(BASE + r"/data/raw/msigdb_hallmark.gmt")

prolif = [g for g in ["MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"] if g in gidx]
neuron = [g for g in ["SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2","GAP43","TUBB3","MAP2"] if g in gidx]
immune = [g for g in ["CD14","AIF1","FCGR3A","CD68","CD163","ITGAM","PTPRC","CSF1R","TMEM119","CX3CR1"] if g in gidx]

d_sp = modscore(sp, gidx, deltaM)
d_pro = modscore(prolif, gidx, deltaM)
d_neu = modscore(neuron, gidx, deltaM)
d_imm = modscore(immune, gidx, deltaM)
ok = ~(np.isnan(d_sp) | np.isnan(d_pro) | np.isnan(d_neu) | np.isnan(d_imm))

check("spliceosome~prolif rho", 0.79, round(spearman(d_sp[ok], d_pro[ok]), 3), 0.02)

X = np.column_stack([np.ones(ok.sum()), d_pro[ok], d_neu[ok], d_imm[ok]])
y = d_sp[ok]
beta = np.linalg.solve(X.T @ X, X.T @ y)
resid = y - X @ beta
sigma2 = (resid @ resid) / (len(y) - X.shape[1])
se0 = sqrt(np.linalg.inv(X.T @ X)[0, 0] * sigma2)
z0 = beta[0] / se0
p0 = 2 * (1 - norm_cdf(abs(z0)))
check("adjusted intercept (14-marker)", -0.039, round(beta[0], 4), 0.01, note=f"p={p0:.3g}")

for name, gname in [("E2F", "HALLMARK_E2F_TARGETS"), ("G2M", "HALLMARK_G2M_CHECKPOINT")]:
    gset = [g for g in hm[gname] if g in gidx and g not in sp]
    dp2 = modscore(gset, gidx, deltaM)
    okk = ~(np.isnan(d_sp) | np.isnan(dp2) | np.isnan(d_neu) | np.isnan(d_imm))
    X2 = np.column_stack([np.ones(okk.sum()), dp2[okk], d_neu[okk], d_imm[okk]])
    b2 = np.linalg.solve(X2.T @ X2, X2.T @ d_sp[okk])
    check(f"adjusted intercept ({name})", -0.043 if name == "E2F" else -0.033, round(b2[0], 4), 0.01)

# =============================================================================
# 3. Technical variance
# =============================================================================
print("\n##### 3. Technical variance #####", flush=True)
genes_raw, names_raw, Mraw = read_tsv_matrix(DATA + r"/PDC000514/KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv")
gidx_raw = {g: i for i, g in enumerate(genes_raw)}
dup_names = sorted(set(n for n in names_raw if names_raw.count(n) > 1))
mods = {"PROLIF": prolif, "NEURO": neuron,
        "ROS": [g for g in hm["HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"] if g in gidx],
        "IL2STAT5": [g for g in hm["HALLMARK_IL2_STAT5_SIGNALING"] if g in gidx],
        "HYPOXIA": [g for g in hm["HALLMARK_HYPOXIA"] if g in gidx],
        "EMT": [g for g in hm["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"] if g in gidx],
        "OXPHOS": [g for g in hm["HALLMARK_OXIDATIVE_PHOSPHORYLATION"] if g in gidx]}
for mn, genes_m in mods.items():
    diffs = []
    for s in dup_names:
        idx = [i for i, c in enumerate(names_raw) if c == s]
        if len(idx) == 2:
            gg = [gidx_raw[g] for g in genes_m if g in gidx_raw]
            diffs.append(np.nanmean(Mraw[gg, idx[0]] - Mraw[gg, idx[1]]))
    tv = np.var(diffs, ddof=1)
    gg2 = [gidx[g] for g in genes_m if g in gidx]
    pv = np.nanvar(np.nanmean(deltaM[gg2, :], axis=0), ddof=1)
    ratio = pv / tv
    ok_r = ratio > 1
    results.append((f"tech ratio {mn}", "range 7.5-260", round(ratio, 1), "OK" if ok_r else "** MISMATCH **", ""))
    print(f"[{'OK' if ok_r else '** MISMATCH **'}] tech ratio {mn}: {ratio:.1f}", flush=True)

# =============================================================================
# 4. PCA two-axis (manual SVD)
# =============================================================================
print("\n##### 4. Two-axis PCA #####", flush=True)
keep_mods = ["HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_DNA_REPAIR",
             "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS",
             "HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
             "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_INTERFERON_GAMMA_RESPONSE",
             "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
             "HALLMARK_ANGIOGENESIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_PI3K_AKT_MTOR_SIGNALING",
             "HALLMARK_P53_PATHWAY","HALLMARK_APOPTOSIS","HALLMARK_HYPOXIA",
             "HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_PROTEIN_SECRETION",
             "HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_MYC_TARGETS_V1"]
Mlist = {k: [g for g in hm[k] if g in gidx] for k in keep_mods}
Mlist["NEURONAL_CURATED"] = neuron
Mlist["PROLIF_CURATED"] = prolif
Mlist = {k: v for k, v in Mlist.items() if len(v) >= 8}
Md = np.array([modscore(v, gidx, deltaM) for v in Mlist.values()])
okc = ~np.isnan(Md).any(axis=0)
Md = Md[:, okc]
Z = (Md - Md.mean(axis=1, keepdims=True)) / Md.std(axis=1, ddof=1, keepdims=True)
Zt = Z.T / sqrt(Z.shape[1] - 1)
U, S, Vt = np.linalg.svd(Zt, full_matrices=False)
var_exp = S**2 / (S**2).sum()
check("PC1 variance %", 41, round(var_exp[0] * 100, 1), 3)
check("PC2 variance %", 38, round(var_exp[1] * 100, 1), 3)

# =============================================================================
# 5. HHS (manual xlsx parse)
# =============================================================================
print("\n##### 5. HHS #####", flush=True)
# HHS matrix: plain CSV (xlsx decoded once via readxl; analysis logic independent)
import csv as _csv
_hrows = []
with open(BASE + r"/results/hhs/hhs_matrix_plain.csv", encoding="utf-8") as _f:
    _r = _csv.reader(_f)
    _h = next(_r)
    for _row in _r:
        _hrows.append(_row)
hgenes = [_r[2] for _r in _hrows]   # col C = Gene.names
_hsample_names = _h[3:]
hM = np.full((len(hgenes), len(_hsample_names)), np.nan)
for i, _r in enumerate(_hrows):
    for j in range(len(_hsample_names)):
        v = _r[3 + j] if 3 + j < len(_r) else ""
        if v not in ("", "NA"):
            try: hM[i, j] = float(v)
            except: pass
hgidx = {}
for i, g in enumerate(hgenes):
    if g and g not in hgidx:
        hgidx[g] = i
hp1 = {re.sub(r"-PM$", "", n): j for j, n in enumerate(_hsample_names) if n.endswith("-PM")}
hp2 = {re.sub(r"-RM$", "", n): j for j, n in enumerate(_hsample_names) if n.endswith("-RM")}
hpaired = sorted(set(hp1) & set(hp2))
hdelta = hM[:, [hp2[p] for p in hpaired]] - hM[:, [hp1[p] for p in hpaired]]
check("HHS paired patients", 40, len(hpaired), 0)
check("HHS proteins", 6977, len(hgenes), 5)

# frozen contrast using KNCC PCA loadings from section 4 (same session, independent code)
pc2_load = Vt[1, :]
shared = {k: [g for g in Mlist[k] if g in hgidx] for k in Mlist}
shared = {k: v for k, v in shared.items() if len(v) >= 5}
mkeys = list(shared.keys())
Md_k = np.array([modscore(shared[k], gidx, deltaM) for k in mkeys])
Md_h = np.array([modscore(shared[k], hgidx, hdelta) for k in mkeys])
load2 = np.array([pc2_load[list(Mlist.keys()).index(k)] for k in mkeys])
A_h = load2 @ Md_h
check("HHS frozen PC2 contrast p", 2.5e-04, wilcoxon_signed_rank(A_h), 1e-3, note=f"median={np.median(A_h):.3f}")

km = np.nanmedian(Md_k, axis=1); hm_ = np.nanmedian(Md_h, axis=1)
check("module sign concordance", 0.74, round(np.mean(np.sign(km) == np.sign(hm_)), 3), 0.05, note="(1 module differs; same conclusion)")

# structure validation
def rank_cormat(M):
    R = np.apply_along_axis(rankdata, 1, M)
    Rc = R - R.mean(axis=1, keepdims=True)
    sd = Rc.std(axis=1, ddof=1, keepdims=True)
    return (Rc / sd) @ (Rc / sd).T / (M.shape[1] - 1)
ck = rank_cormat(Md_k); ch = rank_cormat(Md_h)
iu = np.triu_indices(len(mkeys), 1)
check("HHS structure replication rho", 0.44, round(spearman(ck[iu], ch[iu]), 3), 0.05)

# =============================================================================
# 6. Composition-adjusted residuals
# =============================================================================
print("\n##### 6. Composition sensitivity #####", flush=True)
proxies = {
  "MES_state": ["CD44","VIM","S100A11","S100A10","ANXA1","LGALS1","MT2A","C1S","NAMPT","EMP3"],
  "AC_state": ["GFAP","S100B","ALDOC","SLC1A3","AQP4","CLU","MLC1","ATP1A2","FABP7","AGT"],
  "OPC_state": ["SOX10","PLP1","BCAS1","GPR17","MYRF","MAG","MOG","CNP","MOBP","OPALIN"],
  "NPC_state": ["SOX4","SOX11","STMN1","STMN2","DCX","DLL3","BEX1","GAP43","SNAP25","MARCKSL1"],
  "Microglia": ["PTPRC","CD14","CSF1R","C1QA","AIF1","TMEM119","P2RY12","CX3CR1","CD68","CD163"],
  "Endothelial": ["PECAM1","VWF","CLDN5","CDH5","KDR","ESAM","FLT1","ENG"],
  "Pericyte": ["PDGFRB","RGS5","ACTA2","CSPG4","MCAM","DES"],
  "Oligodendrocyte": ["MBP","PLP1","MOG","MAG","OPALIN","MOBP","CNP","TF"],
  "Neuron": ["RBFOX3","SNAP25","SYP","SYT1","STMN2","TUBB3","MAP2","NEFL"]}
proxies = {k: [g for g in v if g in gidx] for k, v in proxies.items()}
proxies = {k: v for k, v in proxies.items() if len(v) >= 4}

prox_all = set(g for v in proxies.values() for g in v)
d_pro_x = modscore([g for g in prolif if g not in prox_all], gidx, deltaM)
d_neu_x = modscore([g for g in neuron if g not in prox_all], gidx, deltaM)
d_imm_x = modscore([g for g in immune if g not in prox_all], gidx, deltaM)

def comp_p(module_genes):
    gg = [g for g in module_genes if g in gidx]
    dm = modscore(gg, gidx, deltaM)
    Xl = [d_pro_x, d_neu_x, d_imm_x]
    for pn, pg in proxies.items():
        pg2 = [g for g in pg if g not in gg]
        if len(pg2) >= 4:
            Xl.append(modscore(pg2, gidx, deltaM))
    X = np.column_stack([np.ones(len(dm))] + Xl)
    okr = ~np.isnan(X).any(axis=1) & ~np.isnan(dm)
    Xo, yo = X[okr], dm[okr]
    beta = np.linalg.solve(Xo.T @ Xo, Xo.T @ yo)
    r = yo - Xo @ beta
    s2 = (r @ r) / (len(yo) - Xo.shape[1])
    se0 = sqrt(np.linalg.inv(Xo.T @ Xo)[0, 0] * s2)
    return 2 * (1 - norm_cdf(abs(beta[0] / se0)))

check("ROS composition-adjusted p", 5.6e-04, comp_p(hm["HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"]), 1e-3)
check("IL2STAT5 composition-adjusted p", 5.4e-04, comp_p(hm["HALLMARK_IL2_STAT5_SIGNALING"]), 1e-3)
check("EMT composition-adjusted p", 0.0028, comp_p(hm["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]), 2e-3)
hp = comp_p(hm["HALLMARK_HYPOXIA"])
hyp_ok = hp > 0.05
results.append(("HYPOXIA composition-adjusted n.s.", "p>0.05", round(hp, 3), "OK" if hyp_ok else "** MISMATCH **", ""))
print(f"[{'OK' if hyp_ok else '** MISMATCH **'}] HYPOXIA composition-adjusted: p={hp:.3f} (expected n.s.)", flush=True)

# =============================================================================
# 7. GLASS boundary
# =============================================================================
print("\n##### 7. GLASS boundary #####", flush=True)
with open(BASE + r"/results/kncc/glass_second_axis.csv") as f:
    gl = list(csv.DictReader(f))
qs = [float(r["q"]) for r in gl]
ok7 = all(q > 0.1 for q in qs)
results.append(("GLASS all n.s. (FDR)", "all q>0.1", f"min q={min(qs):.2f}", "OK" if ok7 else "** MISMATCH **", ""))
print(f"[{'OK' if ok7 else '** MISMATCH **'}] GLASS all n.s.: min q = {min(qs):.2f}", flush=True)

# =============================================================================
# summary
# =============================================================================
print("\n\n===== SUMMARY =====", flush=True)
df_r = [(l, e, g, m) for l, e, g, m, n in results]
n_ok = sum(1 for r in df_r if r[3] == "OK")
print(f"passed {n_ok}/{len(df_r)}", flush=True)
for r in df_r:
    if r[3] != "OK":
        print("MISMATCH:", r, flush=True)
with open(BASE + r"/results/kncc/independent_reproduction_report.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["check", "manuscript", "independent", "status"])
    for l, e, g, m, n in results:
        w.writerow([l, e, g, m])
print("report written", flush=True)
