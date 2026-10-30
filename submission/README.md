# Computers in Biology and Medicine — 投稿包

生成日期：2026-10-30 | 所有图片经逐张目检，空白页/错误图例/截断坐标轴已修复

## 五类结构

```
submission/
├── Cover_Letter/
│   └── 01_Cover_Letter.docx          # 投稿信（novelty/scope/AI声明）
├── Manuscript/
│   ├── 02_Manuscript.docx            # 主文档（标题页→结论→声明）
│   ├── 03_Declarations.docx          # 声明汇总
│   └── 04_Suggested_Reviewers.docx   # 推荐审稿人
├── Figures/                          # 正文主图（jpg 300dpi / pdf / tif 300dpi）
├── Tables/                           # 正文表格（Table 1 队列、Table 2 靶点分级）
└── Supplementary/
    ├── Supplementary_Figures.pdf     # 补充图合并
    ├── Supplementary_Tables.xlsx     # 25 个补充表
    └── figures/                      # 补充图单独文件（jpg/pdf/tif）
```

## 主图清单（Figures/）与正文对应

| 文件 | 内容 | 对应章节 |
|------|------|---------|
| Figure1_StudyDesign_v16_modified | 研究设计（原图） | Fig 1 |
| fig_copykat_purity | CopyKAT 纯度控制 | Results 2.1 |
| fig_pseudotime | 拟时序轨迹 | Results 2.3 |
| fig_acetyl_spliceosome_heatmap | 剪接体乙酰化位点热图 | Results 2.5 |
| fig_spatial_GBM1~GBM5_2 | 空间 WEBS（6 张，对称色标） | Results 2.7 |
| fig_survival_TCGA_GBM / cgga693 / cgga325 / *IDHwt | WEBS 分组 KM（单页修复版） | Results 2.8 |
| fig_signature_TCGA_GBM / cgga693 / cgga325 | 六基因签名 KM（固定系数） | Results 2.8 |
| fig_signature_lasso | LASSO 路径图 | Results 2.8 |
| fig_depmap_scatter | DepMap 依赖性（MED1/SF3B1/MLLT1/KAT2A） | Results 2.9 |
| fig_gdsc2_scatter | HDAC 抑制剂药敏 | Results 2.9 |
| fig_cmap_drugs | 药物重定位 Top 20 | Results 2.9 |
| fig_cellchat_compare | LLS vs LHS 互作数/强度 | Results 2.10 |
| fig_cellchat_ranknet | 通路信息流排序 | Results 2.10 |
| fig_cellchat_diffnet | 差异互作网络 | Results 2.10 |
| Figure12_WorkingModel_v16 | 工作模型（原图） | Fig 12 |

## 补充图（Supplementary/figures/）

figS_bootstrap / figS_pancancer（5 页）/ figS_performance / figS_sensitivity（修复版）/ figS_shap / figS_specificity

## 作者投稿前需手动补齐

1. **References**：从原稿 BBA_v16_2.docx 复制到 Manuscript 末尾，套 Elsevier numbered 样式
2. **Funding**：02_Manuscript.docx 中填基金号
3. **原论文其余主图**（Figure 2–11 中旧版图）：从 BBA_v16_2.docx / SI_v8_2.docx 导出后放 Figures/
4. **Graphical Abstract**：按约定未生成

## 上传建议（Editorial Manager）

- 图优先上传 tif（300 dpi）；pdf 为矢量备选
- 多页 pdf 转出的 _p1/_p2 文件已消除（修复版均为单页）
- Article Type: Research Article
