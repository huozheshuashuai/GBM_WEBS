# 新增分析结果清单（复现完成之后）

日期：2026-11-05。以下全部为原论文之外的新增分析，用于投稿强化版。

## A. 机制深化
| 分析 | 关键结果 | 文件 |
|------|---------|------|
| 乙酰化位点机制 | 剪接体乙酰化位点 OR=7.36；耦合 ρ=0.62 | results/acetyl_summary.txt |
| 剂量反应 | 43/43 位点随 WEBS 单调递增（Jonckheere） | results/jonckheere_trend.csv |
| 转录-蛋白背离尾部 | p=4.7e-49 | figures_final/fig_discordance_tail.pdf |
| **可变剪接功能验证** | 4,387 个 WEBS 相关 AS 事件（15.6%）；富集于 RNA splicing 本身（剪接级联） | results/splicing_webs_corr.csv, fig_splicing_volcano.pdf, fig_splicing_eventtypes.pdf, fig_splicing_enrichment.pdf |
| **磷酸化蛋白组** | 198/319 剪接体磷酸化位点显著（Fisher OR=3.26, p=3e-25）；SRSF 家族 45/71 显著；剪接激酶 CLK3 ρ=0.49、CDK9 ρ=0.56、CDK12 ρ=0.54 | results/phospho_summary.txt, fig_phospho_volcano.pdf |
| 虚拟 KO（GENIE3 GRN） | 剪接体基因对全部 14 个表观酶 KO 的敏感性 > 其他基因；HDAC2/SIRT2 最强 | results/virtual_ko_summary.txt, fig_virtual_ko.pdf |
| STRING PPI | 121 节点、3349 边、p=1e-64；hub: SNW1/SRSF1/SF3B1 | results/ppi_summary.txt, fig_ppi_network.pdf |

## B. 临床转化
| 分析 | 关键结果 | 文件 |
|------|---------|------|
| 治疗交互 | risk×TMZ p=0.0024；risk×RT p=0.00045（CGGA-693） | results/confounder_subgroup_summary.txt |
| MGMT 联合分层 | 两队列 p<1e-13 | fig_mgmt_joint_cgga325/693.pdf |
| 混杂校正 | CGGA-325 校正后 p=0.0011 仍显著；CGGA-693 校正后减弱 p=0.182（诚实报告） | 同上 |
| Meta 分析 | 合并 HR=1.35, p=0.0028 | fig_meta_forest.pdf |
| **GLASS 复发时间** | 初发 WEBS ρ=+0.327, p=1.1e-06；高 WEBS 中位复发 18 月 vs 低 10 月 | results/glass_recurrence_time.txt, fig_glass_recurrence_time.pdf |

## C. 外部验证
| 队列 | 结果 |
|------|------|
| **GLASS（n=107 原发）** | 六基因签名 **HR=2.14, p=2.0e-06**；KM p=4.2e-05 |
| GLASS 纵向（179 对） | 复发时 WEBS 变化 n.s.（p=0.075） |
| GSE83300（n=50） | HR=2.76 方向一致，p=0.74（功效不足） |

## D. 方法学稳健性
| 分析 | 结果 |
|------|------|
| 6 方法 benchmark | scRNA 拟时序相关 WEBS 最强（0.303）；bulk C-index 非最强（0.489，与边界声明自洽） |
| Bootstrap 稳定性 | HDAC1 100%、LDHB 97%、VEGFA 91% 入选率 |
| 参数敏感性 | maxRank + leave-one-out 稳定 |
| 计算性能 | WEBS vs AUCell 计时/内存报告 |
| 特异性对照 | 100 随机基因组 empirical p=0.010 |

## E. 多组学整合
| 分析 | 结果 |
|------|------|
| SNF 跨组学聚类 | KW p=7.06e-07 |
| WGCNA | ME3 r=0.72（神经谱系模块，诚实报告非剪接体） |
| MOFA2 | Factor 2 ρ=−0.64 |
| 空间自相关 | Moran's I 显著（半数切片） |
| 通路级 SHAP | 完成 |

## F. 泛癌（完成）
- 33 TCGA 癌种 WEBS×剪接体相关 + 签名生存：14/33 显著；GBM 签名 HR=2.88 p=3.2e-07 全泛癌最强：results/pancancer_33_summary.csv

## G. 统计严谨性强化（投稿前）
| 分析 | 关键结果 | 文件 |
|------|---------|------|
| 伪重复控制（pseudobulk+混合模型） | Neftel 梯度患者级 p=7.4e-04 保持；缺氧解耦 p=1 确认；区域/拟时序患者级减弱（诚实降级） | results/pseudobulk_validation.txt, fig_pseudobulk.pdf (S21) |
| PTM 母蛋白丰度校正 | 乙酰化 OR 7.36→2.29 (p=0.003)；磷酸化 OR 3.26→2.95 (p=1.3e-18) | results/ptm_abundance_correction.txt |
| 特异性 10k 表达匹配对照 | p=0.22（推翻原 p=0.010 声明，已诚实改写） | results/specificity_10k_summary.txt, figS_specificity_10k.pdf (S7 替换) |
| Uncoupling 严谨化 | Δρ=0.81 bootstrap CI [0.60,1.03] p<0.0001；LOO 0.788-0.795；控制增殖/免疫 ρ=0.65 | results/uncoupling_rigor.txt, fig_uncoupling_pergene.pdf (S22) |
| WEBS 增量价值 | 控制增殖/MES/缺氧/OPC 后 ρ=0.66 p=5.8e-13 | results/incremental_value.txt |
| 剪接 TME 校正 | 91% 事件保留（3987/4387） | results/splicing_tme_summary.txt |
| 剪接 LGG 跨队列复现 | （进行中） | results/splicing_lgg_summary.txt |
