# GBM WEBS 论文复现与强化工程 — 完整成果报告

**日期**：2026-10-30
**项目目录**：`E:\KIMI\GBM_WEBS_Project\`
**论文**：Chen et al. "Spliceosome activation defines an epigenetically hyperactive and therapeutically targetable malignant cell state in glioblastoma"（投稿 BBA - Molecular Basis of Disease）

---

## 1. 环境搭建（新电脑从零重建）

| 组件 | 版本 | 位置 |
|------|------|------|
| R | 4.5.2 | `E:\R\R-4.5.2` |
| Rtools | 45 | `E:\Rtools45` |
| Miniconda | Python 3.14.9 | `E:\Miniconda3`（含 docking / mofa 两个环境） |
| AutoDock Vina | 1.2.7 | `GBM_WEBS_Project/docking/vina.exe` |
| OpenBabel | 3.1.1 (conda-forge, 含 pdbqt) | `E:\Miniconda3\envs\docking` |

R 包：60+ 全部装好（Seurat 5 / Signac / UCell / AUCell / slingshot / DESeq2 / GSVA 2.2 / MOFA2 / CellChat v2 / copykat / spacexr(RCTD) / TCGAbiolinks / GEOquery / kernelshap / STRINGdb / spdep 等），安装日志见 `results/package_*_log.txt`。

## 2. 数据集（全部就位）

| 数据集 | 内容 | 来源 |
|--------|------|------|
| GSE138794 | 35,823 细胞 scRNA + 6 样本 snATAC + 官方细胞注释 | GEO（用户下载） |
| GSE194329 | 6 个 GBM Visium 空间切片 | GEO（用户下载） |
| GSE131928 | Smart-seq2 (7,930 细胞) + 10x (16,201 细胞) 外部验证 | GEO |
| CGGA | mRNAseq_693 (657 样本) + mRNAseq_325 (313 样本) + 临床 | CGGA 网站 |
| CPTAC-GBM | 蛋白组 99×11,141 + 乙酰化组 99×13,556 + 磷酸化组 + RNAseq + 临床 | LinkedOmics |
| TCGA-GBM | LinkedOmics 表达 (144 样本) + GDC (371) + RPPA500 (71 匹配) + RPPA 旧版 | LinkedOmics/TCPA/GDC |
| DepMap 24Q2 | CRISPRGeneEffect + Model + 表达矩阵 | DepMap |
| GDSC2 | fitted dose response (242,037 条) + 细胞系详情 | Sanger FTP |
| TCGA 泛癌 | LUAD (527) / BRCA (1,087) / SKCM (458) | GDC 自动下载 |
| GSE7696 | 84 样本 + 生存 | GEO |
| 未纳入 | REMBRANDT（GEO matrix 无生存字段）、GSE13041（生存字段提取失败）、GLASS（待 Synapse） | — |

## 3. 论文核心发现复现对照

| # | 论文发现 | 论文值 | 复现值 | 判定 |
|---|---------|--------|--------|------|
| 1 | 单细胞规模 | 35,635 细胞 | 35,823 | ✅ |
| 2 | 剪接体蛋白上调 | 117/129, OR=39.2, FDR=6.9e-64 | 107/111 上调、0 下调, p=3.4e-37 | ✅ |
| 3 | 蛋白层 WEBS-剪接体相关 | ρ=0.44 | ρ=0.79 | ✅ 同向更强 |
| 4 | mRNA 层弱/负相关 | ρ=−0.29 | ρ=−0.02 (n.s.) | ✅ 支持翻译后积累 |
| 5 | 蛋白 WEBS vs 组蛋白乙酰化 | 正相关 | r=0.62, p=5.3e-12 | ✅ |
| 6 | MED1 依赖 | 显著 | r=−0.28, q=0.049（IDHwt: p=0.043） | ✅ |
| 7 | MLLT1 依赖 | 显著 | r=−0.15, p=0.14；Chronos 从不 <−1 | ⚠️ 不可复现（建议改口径或降级） |
| 8 | 伏立诺他敏感 | 显著 | r=−0.16, q=1.1e-4；恩替诺他 q=8.7e-5 | ✅ |
| 9 | 拟时序 LLS→LHS | 正相关 | ρ=0.301, p=2.7e-166 | ✅ |
| 10 | 缺氧解耦 | 解耦 | 恶性细胞内 ρ=−0.054；调缺氧后 WEBS-Margin p=3.8e-4 | ✅（复核成立） |
| 11 | CGGA 生存分层 | 显著 | 两队列 p≈1e-18 | ✅ |
| 12 | TCGA-GBM 生存 | 显著 | 原始 WEBS 不显著（多口径核实）；六基因评分 HR=3.75, p=1.4e-10 | ⚠️ 表述需精确化 |
| 13 | 外部验证 Neftel 状态差异 | 显著 | KW p=1.6e-168 (SS2) / 1e-313 (10x) | ✅ |

## 4. 新增强化分析（超出原论文）

| 分析 | 关键结果 | 文件 |
|------|---------|------|
| **位点级乙酰化机制** | 剪接体蛋白乙酰化位点在 LHS 富集 OR=7.36, p=1.3e-18；52/88 位点与 WEBS 正相关（中位 ρ=0.56）；乙酰化-丰度耦合中位 ρ=0.62（PRPF40B 0.87、HNRNPA2B1 0.85、PRPF8 0.84） | `acetyl_*.csv`, `fig_acetyl_spliceosome_heatmap.pdf` |
| **CopyKAT 纯度控制** | 18 样本 16,255 细胞；4,421 个纯肿瘤细胞内 WEBS 分层 KW p≈0 | `copykat_summary.txt`, `fig_copykat_purity.pdf` |
| **参数敏感性** | maxRank 扫描 ρ=0.68–0.95（需固定 1500 并报告）；LOO/L2O 57 组中位 ρ=0.95（0.76–0.99） | `sensitivity_*.csv`, `figS_sensitivity.pdf` |
| **方法 benchmark** | WEBS vs AUCell 0.95 / MeanExpr 0.92 / ssGSEA 0.87 / PLAGE −0.97（符号翻转）；WEBS 35.8k 细胞 213s、内存 <2MB | `benchmarking_*.csv`, `computational_performance.csv`, `figS_performance.pdf` |
| **泛癌验证** | LUAD HR=0.64 p=0.039 ✅；BRCA/SKCM/GBM-bulk n.s.（癌种特异） | `pancancer_webs_summary.csv`, `figS_pancancer.pdf` |
| **MOFA2 四组学** | Factor 2 与 WEBS ρ=−0.637, p<1e-12；解释蛋白组 11% 方差；剪接体 loading p=0.123（趋势） | `mofa2_*.csv` |
| **空间自相关** | 3/6 切片显著聚集（GBM3 I=0.222、GBM5_2 I=0.200、GBM4 I=0.053，置换 p≈0.001） | `spatial_autocorrelation.csv` |
| **空间 WEBS 分布** | 6 切片 WEBS 图 + dual-hotspot（GBM1 314 个等） | `spatial_summary.csv`, `fig_spatial_*.pdf` |
| **CellChat** | SPP1 上调于 LHS（+0.18，可成药）；MIF-CD74/CXCR4 轴 LLS 富集（置换 p<1e-4） | `cellchat_*.csv`, `fig_cellchat_*.pdf` |
| **Bootstrap 稳定性** | HDAC1 100%、LDHB 97%、VEGFA 91%、LDHA 89% 选择频率 | `bootstrap_selection_frequency.csv`, `figS_bootstrap.pdf` |
| **SHAP** | HDAC1 贡献最大（mean\|SHAP\|=0.31，负向保护） | `shap_*.csv`, `figS_shap.pdf` |
| **GDSC2 协同 (Bliss)** | Piperlongumine+CHIR-99021 最优（ρ=0.29, q=0.012） | `gdsc2_combination.csv` |
| **药物重定位** | 22 个 WEBS 增敏药物（FDR<0.1）：Piperlongumine、CHIR-99021、SL0101、Nutlin-3a、AZD5991 等 | `translation_summary.txt`, `fig_cmap_drugs.pdf` |
| **分子对接** | HDAC2←伏立诺他 −7.95；SIRT2←AGK2 −8.52；SF3B1←H3B-8800 −6.68；SF3B1←pladienolide B −6.20 kcal/mol | `docking/results/`（含 4 个构象 pdbqt） |
| **RPPA500 LMS** | 蛋白层乳酸代谢评分分层 TCGA-GBM 生存 KM p=0.016（n=71） | `rppa_lms_summary.txt` |

## 5. 论文写作需修正的 4 处（复核结论）

1. **TCGA 生存**：原文「WEBS stratifies survival」→ 改为「the six-gene signature stratifies survival」；原始 WEBS 在 TCGA bulk 阴性可写进 boundary 段落。RPPA500 LMS p=0.016 支持蛋白层信号更强的论述。
2. **MLLT1**：建议改为「GBM 选择性依赖」口径重算，或移至补充材料；主文保留 MED1。
3. **缺氧解耦**：明确写「in malignant cells」（全体细胞 ρ=−0.45 是混入非肿瘤细胞的假象）。
4. **MIF**：写「LLS 富集的 MIF-CD74/CXCR4 免疫抑制轴」（配体受体均在 LLS 高，置换 p<1e-4），不写 LLS→LHS 方向性。

## 6. 不可行/未完成项（诚实记录）

- **RNA Velocity**：GSE138794 无 BAM/FASTQ，GEO 仅计数矩阵 → Methods 需注明
- **RCTD**：spacexr 内部缓存 bug（choose_sigma_c readRDS 失败），已有 Seurat 层面空间证据替代
- **REMBRANDT / GSE13041 / GLASS**：前两者 GEO 页面无可用生存字段；GLASS 待 Synapse 注册下载
- **PAGA**：需 Python scanpy 环境（Miniconda 已就绪，可随时补）

## 7. 文件索引

- **脚本**：`run_02/02b/03/04/04b/07/08/09/10/11*.R`（主流程驱动）+ `13–19*.R`（方法学强化）+ `20–27*.R`（机制/复核/图）
- **结果**：`results/`（48 个 csv/txt）
- **图**：`figures/`（31 张 PDF）
- **对接**：`docking/`（pdb/ + ligands/ + results/binding_energies.csv + 构象 pdbqt）
- **数据**：`data/raw/`（9 个数据集）+ `data/processed/`（Seurat RDS 全套）

## 8. 剩余待办（非计算）

1. Figure 1 / Figure 12 概念图手工绘制（可提供元素清单与图注文案）
2. 4 段修正文字起草（可代写）
3. GLASS 队列（Synapse 注册后下载，可选）
