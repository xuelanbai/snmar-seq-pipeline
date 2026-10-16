# DNA TSS 高变 promoter 敏感性分析（步骤 147–150）

## 目的

步骤 122–134 的 TSS±1.5 kb 结果没有显示可解释的 GM12878/PK15
分离。该结果中每个细胞的 WCG promoter 覆盖很低，同时
`max-features-per-context=20000` 使约 6,800 个候选 promoter 几乎全部
进入降维。本阶段只改变预先定义的稀疏度和高变特征过滤，不改变窗口、
物种标签或 UMAP 参数，也不重建 96 孔矩阵。

GM12878/PK15 mapping 标签仅用于结果完成后的评价，不参与细胞过滤、
特征排序、SVD 或 UMAP。

## 固定参数网格

共运行四个作业：

- 全基因组 DNA calls：D20 与 D50；
- 每个 context 的共享 promoter calls：至少 100；
- 每个特征至少在 20 个细胞出现；
- 每个特征至少 100 个 covered calls；
- 选取 Pearson residual variability 最高的 1,000 或 2,000 个特征；
- SVD 30维，UMAP `n_neighbors=30`、`min_dist=0.30`。
- 若高变特征选择后某个细胞在WCG或GCH入选特征中完全无观测，则按
  标签盲法删除该细胞；累计删除上限为该运行候选细胞的1%，超过即报错。

不运行 top 500：在当前稀疏度下容易使低覆盖细胞在入选特征中完全无
观测，并且会扩大参数选择偏差。

## 运行

从项目根目录执行：

```bash
cd /share/org/YZWL/yzwl_wujx/Projects/2607_methylation

bash 03_scripts/05_downstream_analysis/06_dna_tss_umap/148_submit_dna_tss_hvf_4jobs.sh

CONFIRM_SUBMIT_DNA_TSS_HVF=1 \
bash 03_scripts/05_downstream_analysis/06_dna_tss_umap/148_submit_dna_tss_hvf_4jobs.sh
```

步骤 147 是模板，不可直接提交。四个作业全部完成后：

```bash
csub < 03_scripts/05_downstream_analysis/06_dna_tss_umap/150_submit_dna_tss_hvf_summary_4core.csub
```

主要汇总：

```text
08_qc_summary/39_dna_tss_umap/
  dna_tss_hvf_per_run.tsv
  dna_tss_hvf_neighbor_stability.tsv
  dna_tss_hvf_comparison_summary.tsv
```

预先定义的最低支持条件为：任一运行同时满足 WCG silhouette ≥0.10、
WCG 近邻物种一致率高于随机类别比例基线至少0.05。若四个固定运行均
未达到，不继续追逐 UMAP 参数，记录为当前 TSS±1.5 kb 数据下的阴性
结果。

## 扩大 promoter 范围的后续设计

扩大窗口会改变原始计数矩阵，不能复用当前±1.5 kb VDS，也不能覆盖
步骤122–128的结果。高变分析仍无效时，优先建立独立的 TSS±5 kb
结果目录：理论区域长度约为当前的3.3倍，预期提高覆盖，同时仍保留
promoter/TSS解释。TSS±10 kb理论长度约为当前的6.7倍，但更容易合并
相邻基因并混入gene-body或远端调控信号，只作为第二档敏感性对照。

扩大窗口时必须分别使用新索引、新96孔矩阵、新VDS和新结果目录，旧的
±1.5 kb结果保留为基线。

此外，不能把步骤122的“重叠窗口直接合并”原样用于±5/±10 kb：窗口
变宽后，相邻基因的 promoter 更容易合并成多基因 cluster，反而会减少
可配对的单基因共享特征。扩大范围版本将保持每个共享基因一个feature，
重叠区域内的胞嘧啶只分配给最近TSS（距离相同按稳定索引决定），避免
同一call重复进入多个promoter，也避免将邻近基因强制合并。
