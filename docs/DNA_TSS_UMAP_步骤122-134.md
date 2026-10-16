# DNA TSS±1.5 kb 单细胞矩阵与 UMAP（步骤 122–134）

## 目的

这一阶段不替换已经完成的 100 kb 全基因组矩阵，而是增加一套以
TSS 上下游 1.5 kb 为特征的 DNA 分析，用来检查 WCG 甲基化和 GCH
开放性是否能把 GM12878 与 PK15 细胞分开。

主分析不会直接把 human TSS 和 pig TSS 当成两套互斥特征。步骤 122
先找出人、猪中一对一的同名单基因 promoter cluster，再把两物种的信号
投影到同一个共享特征轴。这样可以降低“比对到哪一个参考基因组”本身
造成的物种分离。这里的一对一关系来自大小写归一化的 `gene_name`，是
保守代理，不等同于经过外部数据库确认的严格直系同源关系。

GM12878/PK15 标签来自 DNA 比对至 human/pig 参考序列的 reads 比例，
只在 UMAP 完成后着色和计算近邻一致率，不参与细胞过滤、特征选择、
SVD 或 UMAP；因此它是内部一致性检查，不是独立验证。

## 窗口和缺失值定义

- 每个 gene 的窗口为 `[TSS-1500, TSS+1501)`，包含 TSS 本身以及两侧
  各 1500 bp。
- 到染色体两端时将窗口截断，不生成染色体外坐标。
- 同一 contig 内相互重叠的窗口合并为不重叠 promoter cluster；相邻但
  不重叠的半开区间不合并。
- 单个胞嘧啶最多进入一个 promoter cluster，不因多个重叠转录本重复
  加权。
- WCG 与 GCH 始终分开保存和建模：WCG 表示内源甲基化，GCH 表示
  GpC 酶标记的开放性。
- 矩阵中缺失的 cell-feature 项表示无覆盖（NA），降维时编码为该特征
  期望值对应的 Pearson residual 0，不解释为 0% 甲基化或 0% 开放性。

## 全量运行顺序

以下命令都从项目根目录运行：

```bash
cd /share/org/YZWL/yzwl_wujx/Projects/2607_methylation
```
### 1. 生成 TSS/promoter 索引

```bash
csub < 03_scripts/05_downstream_analysis/06_dna_tss_umap/124_submit_dna_tss_promoter_index_4core.csub
```

主要输出：

```text
10_downstream_analysis/08_dna_tss_promoters/
  dna_tss_promoter_index.tsv.gz
  dna_tss_conservative_ortholog_index.tsv.gz

08_qc_summary/38_dna_tss_promoters/
  dna_tss_promoter_index_summary.tsv
```

确认 summary 中 `status=validated` 后进入下一步。

### 2. 96 个物理孔并行建立 TSS 稀疏矩阵

先预览，不提交：

```bash
bash 03_scripts/05_downstream_analysis/06_dna_tss_umap/126_submit_dna_tss_full_96jobs.sh
```

确认预览后提交 96 个独立作业：

```bash
CONFIRM_SUBMIT_DNA_TSS=1 \
bash 03_scripts/05_downstream_analysis/06_dna_tss_umap/126_submit_dna_tss_full_96jobs.sh
```

每个作业申请 8 CPU。脚本支持断点续跑；已有当前版本的 validated
结果会跳过。步骤 125 是带占位符的模板，不可直接 `csub`。

查看完成数：

```bash
find 08_qc_summary/38_dna_tss_promoters/well_matrices \
  -maxdepth 1 -type f -name '*.dna_tss_promoters.summary.tsv' | wc -l
```

### 3. 合并 96 孔并建立完整 HDF5 VDS

96 个孔均完成后：

```bash
csub < 03_scripts/05_downstream_analysis/06_dna_tss_umap/128_submit_dna_tss_complete_merge_16core.csub
```

主要输出：

```text
10_downstream_analysis/08_dna_tss_promoters/
  dna_tss_complete_matrix.vds.h5

08_qc_summary/38_dna_tss_promoters/
  dna_tss_complete_summary.tsv
  dna_tss_well_manifest.tsv
```

VDS 依赖 96 个孔的源 HDF5 文件，后续不要单独移动或删除
`well_matrices/`。

### 4. 并行运行 D20 与 D50 两套 TSS UMAP

先预览：

```bash
bash 03_scripts/05_downstream_analysis/06_dna_tss_umap/132_submit_dna_tss_umap_2jobs.sh
```

确认后提交两个 32 CPU 作业：

```bash
CONFIRM_SUBMIT_DNA_TSS_UMAP=1 \
bash 03_scripts/05_downstream_analysis/06_dna_tss_umap/132_submit_dna_tss_umap_2jobs.sh
```

两套分析均要求共享 promoter 特征中的 WCG 和 GCH 覆盖各不少于 25
calls；两者仅全基因组 DNA 深度阈值不同：

- D20：WCG+GCH 全基因组 covered calls ≥20,000；
- D50：WCG+GCH 全基因组 covered calls ≥50,000。

主要输出：

```text
10_downstream_analysis/09_dna_tss_umap/
  ortholog_single_gene_D20k_TSSctx25/
  ortholog_single_gene_D50k_TSSctx25/

08_qc_summary/39_dna_tss_umap/
  ortholog_single_gene_D20k_TSSctx25/dna_tss_umap_summary.tsv
  ortholog_single_gene_D50k_TSSctx25/dna_tss_umap_summary.tsv
```

每套结果包含 WCG、GCH 和 joint 三个 UMAP、坐标表、入选特征表及
可重复使用的 SVD/UMAP `npz`。

### 5. 汇总并比较 D20/D50 稳定性

两套 UMAP 均为 validated 后：

```bash
csub < 03_scripts/05_downstream_analysis/06_dna_tss_umap/134_submit_dna_tss_umap_summary_4core.csub
```

主要输出：

```text
08_qc_summary/39_dna_tss_umap/
  dna_tss_umap_per_run.tsv
  dna_tss_umap_D20_D50_stability.tsv
  dna_tss_umap_comparison_summary.tsv
```

稳定性表在 D20/D50 共有细胞中比较 15 近邻集合的 overlap 和 Jaccard。
最终判断不能只看 UMAP 图形是否出现两个团，还要同时查看：

1. `cell_line_knn_agreement` 是否明显高于类别比例造成的随机基线；
2. `cell_line_silhouette` 是否为稳定的正值；
3. D20 与 D50 的近邻结构是否稳定；
4. 分离是否由少量异常孔、覆盖深度或 mapping confidence 驱动；
5. WCG、GCH 和 joint 是否给出相互支持的结果。

## 代码编号

```text
122  准备 TSS/promoter 与共享 gene-name 特征索引
123  单孔 TSS 稀疏矩阵构建程序
124  提交步骤 122
125  单孔作业模板（不可直接提交）
126  预览/提交 96 个单孔作业
127  合并并验证完整 TSS VDS
128  提交步骤 127
129  构建 WCG、GCH、joint TSS UMAP
130  D20 UMAP 作业
131  D50 UMAP 作业
132  预览/提交 D20 与 D50 两个作业
133  汇总 D20/D50 结果和近邻稳定性
134  提交步骤 133
```
