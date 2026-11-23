# Kim et al. 2024 (Cancer Cell) 创新性审计

**审计对象**: Kim KH et al. "Integrated proteogenomic characterization of glioblastoma evolution." Cancer Cell 2024 (PMC10939876; 数据: PDC000514/515, PRJNA1051047)
**审计方法**: 全文 Results 逐节 + 8 张主图图注 + 关键词系统扫描（splice/splicing/acetyl*/ubiquitin/KSEA/WGCNA/MOFA/heterogen*/predict/patient-specific 等）+ GitHub 代码库（ablancsong/Longitudinal_GBM, lucgar/MAKINA, lucgar/GBMstates）
**审计日期**: 2026-11-16

---

## 一、原研究已回答的问题（划掉，不能作为我们的创新点）

| # | 内容 | 原文位置 |
|---|---|---|
| 1 | 纵向基因组演化轨迹（突变稳定性/丢失/新生） | Fig 1 |
| 2 | 复发期超突变表型（8/122 患者）+ MMR 基因改变 + TMZ 突变签名 + 新抗原 | Fig 1, S1 |
| 3 | EGFR 基因组改变与复发时 EGFR 蛋白/磷酸化下降（唯一提及"翻译后机制"处） | Fig 2 |
| 4 | 蛋白组学显示复发期神经元活性上升 | Fig 3 |
| 5 | 多组学亚型分类（PPR/NEU/GPM/MTC）+ 复发期 NEU 亚型富集 + 亚型转换频率 | Fig 4 |
| 6 | 无监督多组学聚类（C1/C2/C3 三簇） | Fig 4 |
| 7 | C2 亚型复发时 MGMT 蛋白升高 + MMR 蛋白下降 → TMZ 耐药（蛋白层面失活机制） | Fig 4, S5 |
| 8 | MGMT 甲基化 vs MGMT 蛋白对生存的不同关联 | Fig S5 |
| 9 | 单细胞（4 对配对）：神经元状态为肿瘤细胞内在 + TME 中少突/神经元增加、髓系减少 | Fig 5 |
| 10 | 突触活性的功能实验验证（湿实验） | Fig 6 |
| 11 | 神经元转变的体内模型（PDC + TMZ 处理纵向采样） | Fig 7 |
| 12 | BRAF 作为治疗靶点 + vemurafenib PDX 前临床验证 | Fig 8 |
| 13 | RTK-RAS/PI3K/p53/细胞周期等核心致癌通路的纵向稳定性 | Fig 1-2 |

**注意**: 原文 Figure 2 中"reduction of EGFR abundance and phosphorylation ... suggesting that post-translational mechanisms operate during evolution"是全文最接近我们主题的一句话，但仅限 EGFR 个案，未展开。

## 二、原研究未回答的问题（关键词零命中确认）

| # | 未做内容 | 证据 | 我们的可行性 |
|---|---|---|---|
| A | **剪接体/可变剪接的任何分析** | "splice"/"splicing" 在 Results 零命中 | 已做（但 H1 被否定，只能当案例） |
| B | **患者级演化异质性/个体化轨迹分析** | heterogen*/patient-specific/inter-tumor/variance 零命中 | ✅ 数据支持（105 对） |
| C | **T1 分子状态预测 T2 演化终点** | "predict" 仅出现在 inferCNV 上下文 | ✅ 数据支持 |
| D | **系统性激酶活性推断（KSEA/激酶-底物网络）** | KSEA 零命中；kinase 仅在通路名称中出现 | ✅ 数据支持（89 对磷酸化），需底物库 |
| E | **蛋白丰度 vs 磷酸化信号的联合演化分解** | 仅 EGFR 个案 | ✅ 数据支持 |
| F | **乙酰化/泛素化等 PTM 层面分析** | acetylat*/ubiquitin 零命中 | ⚠️ KNCC 无乙酰化组；CPTAC 有（横断面） |
| G | **蛋白-mRNA 配对失配的系统量化** | 未做（WTS 与蛋白组在同批患者但未做失配分析） | ⚠️ 需 PRJNA1051047 RNA（991GB 未下载） |
| H | **竞争假说裁决框架（细胞状态 vs 独立调控）** | 未做（原文未区分状态耦合与独立效应） | ✅ 已建立（脚本 91/93/94） |
| I | 生存结局与分子演化轨迹的关联 | 仅 MGMT 相关生存分析 | ✅ 部分已做（需 landmark 严谨化） |
| J | 空间层面的复发验证 | 无空间组学 | ⚠️ GSE254875 GeoMx 可用（147 ROI，需核对结构） |

## 三、审计结论

1. **原研究的核心叙事是"神经元样转变 + BRAF 靶点"**，所有纵向蛋白组/磷酸化组分析都服务于这条线。
2. **我们的竞争假说裁决框架（B→H）在原研究中完全不存在**——这是真正的方法学增量。
3. **患者异质性（B/C）+ 基线预测演化（C）+ 蛋白-磷酸化联合演化（D/E）** 是三个确认空白且数据可及的方向。
4. 剪接体方向已被我们自己的数据否定，只能作为竞争假说框架的演示案例。
5. 原研究有湿实验（Fig 6-8），我们纯计算——定位必须避开与其实验发现直接竞争。

## 四、对新论文设计的含义

可行的差异化定位：**"GBM 复发演化的患者异质性：哪些患者沿哪条路径演化，初发状态能否预测"** —— 原研究回答"平均而言复发时发生什么"，我们回答"不同患者的演化如何不同、能否预测、哪些变化独立于细胞状态转换"。
