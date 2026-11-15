# 03 分析计划（冻结版）

**项目**: GBM Recurrence Evolution Project
**冻结日期**: 2026-11-16
**原则**: 先验证核心现象，再做预测和扩展；不下载大规模原始 RNA；不预设亚型；允许阴性结果终止模块

---

## 一、冻结的核心模块（7 个，不再增加）

| 模块 | 定义 | 角色 |
|---|---|---|
| PROLIF | 14-marker 增殖核心集（MKI67/PCNA/TOP2A/MCM2-7/AURKB/CCNB1/CDK1/BUB1/TTK） | 主要状态轴 |
| NEURO | 10-marker 神经元集（SYP/SNAP25/SYT1/NEFL/NEFM/RBFOX3/STMN2/GAP43/TUBB3/MAP2） | 主要状态轴 |
| ROS | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | B 类候选（独立上调） |
| IL2STAT5 | HALLMARK_IL2_STAT5_SIGNALING | B 类候选 |
| HYPOXIA | HALLMARK_HYPOXIA | PC2 轴组成 |
| EMT | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 状态/组成轴 |
| OXPHOS | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 代谢轴 |

敏感性分析允许使用 Hallmark E2F/G2M 作为增殖替代定义。

## 二、冻结的主要模型

1. **配对变化检验**：每模块 per-patient delta（T2−T1），配对 Wilcoxon，BH-FDR
2. **竞争假说模型**：ΔM ~ ΔPROLIF* + ΔNEURO + ΔIMMUNE（不中心化；截距 = 状态独立效应；预测变量剔除目标模块基因；IMMUNE = 10-marker 集 CD14/AIF1/FCGR3A/CD68/CD163/ITGAM/PTPRC/CSF1R/TMEM119/CX3CR1）
3. **异质性报告**：每模块报告 平均Δ + 95%CI、患者间 SD、方向一致率、患者 bootstrap（2000×）、技术方差上界比较（重复进样）
4. **多重检验**：BH-FDR，按模块族分别校正
5. **排除标准**：主分析 = 全部 105 对；敏感性 = 同 plex 100 对 / 排除 2 名儿科患者（GBM0065, GBM0116）/ 机构内 / leave-one-institution-out / 缺失率阈值变体

## 三、冻结的分析人群

| 人群 | n | 用途 |
|---|---|---|
| 蛋白组配对 | 105 | 主分析 |
| 磷酸化配对 | 89 | 联合演化 |
| RNA 配对（ENA 公开） | 76 | 条件性 RNA 模块 |
| 三层齐全 | 70 | 条件性三层分析 |
| 同 plex 蛋白配对 | 100 | 批次敏感性主对照 |
| 排除儿科后蛋白配对 | ~103 | 敏感性 |

## 四、连续模型优先

患者演化以连续模块评分和 PCA 轴描述。**不预设亚型数量**；只有当离散分群在 bootstrap + 不同特征集 + 不同模块定义下均稳定时才考虑提出亚型（当前 silhouette ≤0.32，暂不支持）。

## 五、条件性模块（允许 No-Go 终止）

| 模块 | 启动条件 | No-Go 处理 |
|---|---|---|
| 基线预测（M0/M1/M2 嵌套 CV） | 核心异质性与联合演化稳定后 | 已完成试点：M2 无增益（perm p 0.42-0.87）→ **已 No-Go，不进入论文主图** |
| RNA-蛋白配对 | 作者提供处理后基因矩阵 | 删除模块 |
| 剪接输出 | 有可靠 junction/PSI 数据 | 不做 de novo 分析 |
| 空间（GSE254875 GeoMx） | 核对患者/ROI/区域结构后可行 | 不纳入核心证据链 |
| PDX | 确认公开时间序列组学存在 | 不纳入 |
| KSEA/激酶 | 校正后位点有信号才做深入 | **已 No-Go（校正后 0/74）** |

## 六、论文产出门禁

- 六张主图不预设；由通过 Go 的模块决定
- 阴性结果（预测、磷酸化、DepMap/GDSC2）如纳入则以"边界结果"形式进入 Discussion
- CPTAC 仅作横向支持，不冒充纵向验证
- 期刊目标由最终核心发现强度决定，不预设 15 分
