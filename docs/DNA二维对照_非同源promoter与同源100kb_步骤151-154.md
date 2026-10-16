# DNA 甲基化 UMAP 的二维对照设计

## 目的

用同一批细胞、相同深度门槛和相同 WCG 指标，区分两类问题：

1. 区域定义：100 kb bin 或 TSS promoter（±1.5 kb）。
2. 跨物种坐标：物种独立坐标或人猪共享坐标。

完整设计为：

| 区域 | 物种独立坐标（reference-aware） | 人猪共享坐标 |
|---|---|---|
| 100 kb | 已有 human bins + pig bins 结果 | 步骤 135–146，同源/共线 100 kb |
| promoter | 步骤 151–154，human promoters + pig promoters | 已有步骤 129–134，单基因同名 promoter 代理轴 |

主结论以 WCG 为准。GCH 和 joint 仅作为开放性及联合信号的补充分析。

## 步骤 151–154：物种独立 promoter

该流程复用步骤 122–128 已生成的完整 promoter 矩阵，不重复扫描 96 孔。
human 和 pig promoter 保持为不同特征轴，最多分别选择 10,000 个 WCG
特征和分别 10,000 个 GCH 特征；这一设计与旧的 human bins + pig bins
100 kb UMAP 对应。

预览两个作业：

```bash
cd /share/org/YZWL/yzwl_wujx/Projects/2607_methylation
bash 03_scripts/05_downstream_analysis/06_dna_tss_umap/152_submit_dna_tss_reference_aware_2jobs.sh
```

提交 D20 和 D50：

```bash
CONFIRM_SUBMIT_DNA_TSS_REFERENCE_AWARE=1 \
bash 03_scripts/05_downstream_analysis/06_dna_tss_umap/152_submit_dna_tss_reference_aware_2jobs.sh
```

两个作业验证完成后汇总：

```bash
csub < 03_scripts/05_downstream_analysis/06_dna_tss_umap/154_submit_dna_tss_reference_aware_summary_4core.csub
```

结果目录：

```text
10_downstream_analysis/09_dna_tss_umap/
  promoter_cluster_D20k_TSSctx25_reference_aware/
  promoter_cluster_D50k_TSSctx25_reference_aware/

08_qc_summary/39_dna_tss_umap/reference_aware_audit/
  dna_tss_reference_aware_audit_summary.tsv
  dna_tss_reference_aware_per_run_embedding.tsv
```

汇总会同时计算全体细胞及仅后 48 个混合孔内的 GM12878/PK15 分离度。
后者是主审计范围，因为前 48 孔只有 GM12878，孔组与物种标签存在混杂。

## 步骤 135–146：同源 100 kb

按 `docs/DNA同源100kb_Seurat试验_步骤135-146.md` 运行。主比较采用：

- D20 与 D50 两套深度门槛；
- `pearson_residual` 作为与现有 Python UMAP 最接近的主矩阵；
- `ratio_paper` 和 `beta_posterior` 作为敏感性分析；
- `coverage_control` 仅用于检验覆盖度是否足以制造类似分群，不能解释为甲基化结果。

## 解释边界

- 物种独立坐标出现分离，可能同时来自甲基化、参考序列身份和结构性缺失。
- 共享 promoter 或同源 100 kb 出现分离，才更接近“相同功能/基因组轴上的甲基化差异”。
- mapping 标签只在 UMAP 完成后着色和评价，不进入矩阵、降维或邻居图。
- 必须报告仅混合孔内的结果，避免把纯 GM12878 孔与混合孔的批次差异误认为物种分离。
