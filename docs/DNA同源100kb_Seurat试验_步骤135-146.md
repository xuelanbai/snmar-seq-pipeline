# DNA 同源 100 kb WCG–Seurat 试验（步骤 135–146）

## 目的

这套试验不替换原有 reference-aware 100 kb UMAP，而是增加一条更严格的
“甲基化特征聚类”路线：human 与 pig 不再作为两套互斥特征，而是通过
GRCh38–Sscrofa11.1 reciprocal-best chain 投影到同一个 orthobin 轴。细胞的
GM12878/PK15 标签只在 Seurat 完成聚类和 UMAP 后加入，用于事后验证。

主分析只使用 WCG 内源甲基化。GCH 开放性和歧义 GCG 不进入本轮聚类。

## 共享特征定义

1. 将 human GRCh38 按固定坐标划分为非重叠 100 kb bin；终端不足100 kb的
   bin可以进入候选，但仍需通过其真实宽度对应的覆盖比例检查。
2. 使用 UCSC `hg38.susScr11.rbest.chain.gz` 将每个human bin投影到pig。
3. 要求单条reciprocal-best chain覆盖至少80%的human bin，pig跨度为human
   bin宽度的0.70–1.30倍。
4. 发生pig区间重叠时，仅保留比对比例/比对碱基数更高的一条，保证一个
   cytosine不会进入两个orthobin。
5. 分别读取human和pig参考基因组，统计两个物种对应区间的可分析WCG：
   plus `[AT]CG` 加上minus `CG[AT]`。`GCG`不计入。
6. human和pig两侧参考WCG都不少于20个时才保留该orthobin。

这里的“参考WCG≥20”是结构过滤。cell-bin甲基化比例的分母仍是该细胞在
该orthobin内实际覆盖的WCG calls，而不是参考基因组WCG总数。

## 四套同轴矩阵

所有矩阵具有完全相同的cells和orthobins：

| 矩阵 | 值 | 无/低覆盖处理 | 用途 |
|---|---|---|---|
| `ratio_paper` | converted/covered WCG | `<5 calls`时使用该orthobin pooled mean | 复现NBT/Seurat思路，主比较 |
| `beta_posterior` | beta-binomial posterior / cell全局WCG均值 | 无覆盖为中性值1 | 稀疏数据稳健主分析 |
| `pearson_residual` | 覆盖校正残差 | 无覆盖为期望值0 | 第二种稳健表示 |
| `coverage_control` | log1p(covered WCG calls) | 无覆盖为0 | 负对照，判断聚类是否由覆盖结构驱动 |

`ratio_paper`在Seurat中执行`NormalizeData(scale.factor=10000)`；其余三套
已经是处理后的数据，直接`ScaleData → PCA → FindNeighbors →
FindClusters → RunUMAP`。

## 服务器运行顺序

都从项目根目录开始：

```bash
cd /share/org/YZWL/yzwl_wujx/Projects/2607_methylation
```

### 0. 将新代码同步到服务器

目录应为：

```text
03_scripts/05_downstream_analysis/07_dna_ortholog_100kb/
```

本阶段尚未验证前不推送 GitHub；可以先通过当前文件传输方式同步到项目。

### 1. 建立human–pig同源100 kb索引

官方chain约79 MB。步骤138可自动下载；下载文件保存在项目参考目录：

```bash
AUTO_DOWNLOAD_CHAIN=1 \
csub < 03_scripts/05_downstream_analysis/07_dna_ortholog_100kb/138_submit_dna_ortholog_100kb_index_8core.csub
```

如果计算节点不能访问外网，在登录节点手动下载：

```bash
mkdir -p 02_reference/human_pig_synteny
curl -L --fail -C - \
  -o 02_reference/human_pig_synteny/hg38.susScr11.rbest.chain.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsSusScr11/reciprocalBest/hg38.susScr11.rbest.chain.gz
```

主要检查：

```text
08_qc_summary/40_dna_ortholog_100kb/
  dna_ortholog_100kb_index_summary.tsv

10_downstream_analysis/10_dna_ortholog_100kb/
  dna_ortholog_100kb_physical_features.tsv.gz
  dna_ortholog_100kb_axis.tsv.gz
```

必须确认summary中：

```text
status = validated
orthobin_axes > 0
physical_features = 2 × orthobin_axes
```

同时记录：低比对比例、异常pig跨度、pig区间重叠和reference-WCG过滤分别
去掉了多少bin。若orthobin过少，先检查contig别名与FASTA是否完全匹配，
不要立即放宽阈值。

### 2. 96孔并行重汇总单碱基WCG/GCH

先预览：

```bash
bash 03_scripts/05_downstream_analysis/07_dna_ortholog_100kb/139_submit_dna_ortholog_100kb_full_96jobs.sh
```

确认后提交：

```bash
CONFIRM_SUBMIT_DNA_ORTHO100=1 \
bash 03_scripts/05_downstream_analysis/07_dna_ortholog_100kb/139_submit_dna_ortholog_100kb_full_96jobs.sh
```

步骤137是模板，不可直接`csub`。96个作业完成后检查：

```bash
find 08_qc_summary/40_dna_ortholog_100kb/well_matrices \
  -maxdepth 1 -type f -name '*.dna_ortholog_100kb.summary.tsv' | wc -l
```

应为96。

### 3. 并行导出D20和D50的四套Seurat矩阵

先预览：

```bash
bash 03_scripts/05_downstream_analysis/07_dna_ortholog_100kb/143_submit_dna_ortholog_100kb_matrix_2jobs.sh
```

确认后：

```bash
CONFIRM_SUBMIT_DNA_ORTHO100_MATRIX=1 \
bash 03_scripts/05_downstream_analysis/07_dna_ortholog_100kb/143_submit_dna_ortholog_100kb_matrix_2jobs.sh
```

两个run：

```text
Dtotal20k_orthoWCG100_call5_top2000
Dtotal50k_orthoWCG100_call5_top2000
```

每个run输出：

```text
10_downstream_analysis/11_dna_ortholog_100kb_seurat_input/<run>/
  ratio_paper.features_by_cells.mtx.gz
  beta_posterior.features_by_cells.mtx.gz
  pearson_residual.features_by_cells.mtx.gz
  coverage_control.features_by_cells.mtx.gz
  features.tsv.gz
  cells.tsv.gz

08_qc_summary/41_dna_ortholog_100kb_seurat_input/<run>/
  dna_ortholog_100kb_seurat_input_summary.tsv
```

候选orthobin还必须在human参考侧和pig参考侧分别至少出现在2个cell、分别
累计至少20个WCG covered calls；这是无标签的参考侧支持检查，防止只有一侧
实际有数据的orthobin进入高变特征。

先检查retained cells、eligible orthobins和selected variable orthobins，确认
矩阵维度不是0且D50是D20的子集，再运行Seurat。

### 4. 确认R/Seurat环境

先在登录节点检查当前`snTAPS_analysis`：

```bash
source /share/apps/anaconda3/etc/profile.d/conda.sh
conda activate snTAPS_analysis
command -v Rscript
Rscript -e 'packageVersion("Seurat"); packageVersion("Matrix")'
```

如果Seurat实际安装在另一个既有环境，提交时设置：

```bash
export DNA_SEURAT_ENV=实际环境名
```

本阶段不要求创建新环境。

### 5. 提交2个深度×4种矩阵，共8个Seurat作业

先预览：

```bash
bash 03_scripts/05_downstream_analysis/07_dna_ortholog_100kb/146_submit_dna_ortholog_100kb_seurat_8jobs.sh
```

确认后：

```bash
CONFIRM_SUBMIT_DNA_ORTHO100_SEURAT=1 \
bash 03_scripts/05_downstream_analysis/07_dna_ortholog_100kb/146_submit_dna_ortholog_100kb_seurat_8jobs.sh
```

初次运行不回归深度。若结果显示methylation UMAP仍与WCG calls梯度高度重合，
再额外提交深度回归敏感性分析：

```bash
REGRESS_DEPTH=1 CONFIRM_SUBMIT_DNA_ORTHO100_SEURAT=1 \
bash 03_scripts/05_downstream_analysis/07_dna_ortholog_100kb/146_submit_dna_ortholog_100kb_seurat_8jobs.sh
```

## 判读原则

不能以“图上出现两团”作为唯一标准。优先比较：

1. `ratio_paper`、`beta_posterior`和`pearson_residual`是否给出相似邻域；
2. GM12878/PK15 mapping标签是否在完成聚类后富集到不同cluster；
3. `coverage_control`是否也产生几乎相同的分离；若相同，主要信号可能是
   物种相关覆盖差异，不足以支持甲基化特征分离；
4. D20和D50中共同细胞的结果是否稳定；
5. 结果是否由少数物理孔或WCG深度梯度驱动；
6. 最后再用RNA严格/放宽标签做独立的post-hoc验证。

只有methylation矩阵稳定分离，而coverage-only对照不能复现，才支持
“仅依靠同源区域甲基化特征区分GM12878和PK15”。

## 编号

```text
135  建立reciprocal-best同源100 kb轴并统计两侧参考WCG
136  单孔单碱基calls重新汇总到同源物理区间
137  单孔作业模板（不可直接提交）
138  提交同源轴/index作业
139  预览/提交96个单孔作业
140  构建四套同轴Seurat输入矩阵
141  D20矩阵作业
142  D50矩阵作业
143  预览/提交D20与D50矩阵作业
144  单个矩阵的Seurat聚类与UMAP
145  Seurat作业模板（不可直接提交）
146  预览/提交2×4共8个Seurat作业
```
