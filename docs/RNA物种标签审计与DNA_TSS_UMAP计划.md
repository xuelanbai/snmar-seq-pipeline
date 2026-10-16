# RNA物种标签审计与DNA TSS UMAP计划

## 当前问题

RNA前48个物理孔按实验设计为GM12878-only。现有汇总中，GeneFull分配到人基因的reads占优，但按`(物理孔, CB, gene, exact UMI)`折叠以后，pig exact UMI和pig细胞标签反而占优。计数代码未发现human/pig字段互换，因此在解释RNA UMAP或把RNA标签连接到DNA之前，必须先定位差异发生在哪些细胞、哪些孔以及哪一个计数层。

`UMI>=500且gene>=200`仅保留3个RNA细胞，不适合降维。当前诊断阈值采用`UMI>=25且gene>=10`，预计保留约4,009个细胞；这只是诊断和首轮RNA UMAP候选阈值，最终阈值仍以逐细胞分布和结果稳定性为准。

## 步骤111–112：RNA read–UMI物种标签审计

主要输入：

```text
08_qc_summary/16_rna_reverse_genefull_full/rna_reverse_cell_qc.tsv.gz
```

脚本分别保留两种标签，不将它们混成一个标签：

- read层：`human_assigned_reads`与`pig_assigned_reads`；
- UMI层：`human_exact_umis`与`pig_exact_umis`；
- 两层均要求主要物种比例至少90%，否则标为`ambiguous`；
- 对高置信标签记录`concordant`或`discordant`；
- 分别计算human和pig的`assigned reads / exact UMI`，检查是否存在明显不对称折叠。

提交命令：

```bash
cd /share/org/YZWL/yzwl_wujx/Projects/2607_methylation
csub < 03_scripts/05_downstream_analysis/05_rna_umap/112_submit_rna_read_umi_cellline_audit_8core.csub
```

主要输出：

```text
10_downstream_analysis/04_rna_cell_line_audit/
  rna_read_umi_cellline_labels.tsv.gz
  rna_read_umi_cellline_audit.png
  rna_read_umi_cellline_audit.pdf

08_qc_summary/34_rna_cell_line_audit/
  rna_read_umi_cellline_audit_summary.tsv
  rna_read_umi_cellline_audit_per_group.tsv
  rna_read_umi_cellline_audit_per_well.tsv
```

只有当GM12878-only孔的read层和UMI层都符合预期且两层冲突率可接受时，RNA物种标签才进入DNA–RNA连接。否则RNA表达UMAP仍可生成，但物种着色必须同时展示read标签、UMI标签和物理孔设计，不能把任一可疑标签当作独立真值。

## RNA UMAP的下一步

审计通过后从完整raw exact-UMI矩阵选细胞，进行library-size归一化、`log1p`、高变基因选择、PCA和UMAP。首轮至少并排显示：

1. 物理孔设计（GM12878-only或GM12878/PK15 mix）；
2. read-derived物种标签；
3. exact-UMI-derived物种标签；
4. UMI数和gene数。

RNA UMAP本身不按物种预先拆分；物种标签仅在嵌入完成后着色。

## DNA TSS上下游1.5 kb UMAP

这条路线值得做，但要避免将“比对到human还是pig参考序列”误当成“甲基化模式分群”。如果直接把human和pig各自的TSS作为互斥特征，缺失模式本身就足以区分物种，不能证明甲基化特征有区分能力。

建议分两层输出：

1. 常规TSS矩阵：按参考基因组定义TSS上下游1.5 kb，WCG为主分析，GCH为开放性对照；窗口在染色体边界处截断，不生成染色体外坐标。
2. 严格跨物种比较：使用human–pig一对一直系同源基因，将两物种各自TSS信号投影到同一“ortholog promoter”特征轴，再做UMAP。这样更接近“依靠甲基化特征区分GM12878和PK15”，而不是依靠参考基因组覆盖差异。

同一物种内相互重叠的TSS窗口需要明确去重策略。用于聚类时优先构建promoter cluster，避免同一个WCG/GCH位点因多个相邻转录本被重复加权；同时保留成员gene/transcript列表，便于解释。

最终验证顺序为：DNA特征独立生成UMAP，之后才通过`物理孔|CB`连接已经审计通过的RNA标签并着色。RNA标签不参与DNA特征选择、PCA或UMAP拟合。

## 步骤111结果与后续严格RNA证据重计数

步骤111在`UMI>=25且gene>=10`下保留4,009个细胞。前48个GM12878-only孔中的结果为：

- human assigned-read加权比例：88.51%；
- human exact-UMI加权比例：23.37%；
- human中位assigned reads/exact UMI：28.10；
- pig中位assigned reads/exact UMI：1.79；
- read层human高置信细胞比例：29.24%；
- UMI层human高置信细胞比例：3.77%。

所有细胞均满足`human_exact_umis <= human_assigned_reads`和`pig_exact_umis <= pig_assigned_reads`，所以不是UMI数超过来源reads的计数越界。异常来自human reads被强烈折叠，而低重复pig背景在exact-UMI层占据较大权重。当前read标签和UMI标签均不能作为DNA UMAP的独立真值。

此前BAM审计还显示GM12878-only孔中仅约0.234%的eligible human reads和0.113%的eligible pig reads带剪接`N` CIGAR。结合RNA严格拆分率偏低及已识别DNA结构残留，不能直接使用现有GeneFull矩阵进行RNA UMAP。

步骤113–117不覆盖旧矩阵，而是从96个现有RNA BAM各扫描一次，同时重算：

1. `genefull`：复现旧GeneFull口径，用于数目守恒；
2. `exon_reverse`：reverse-strand且只命中一个基因exon的reads/UMI；
3. `junction_genefull`：GeneFull分配成功且CIGAR包含`N`的reads/UMI。

先预览4个32核作业：

```bash
bash 03_scripts/05_downstream_analysis/05_rna_umap/115_submit_rna_expression_evidence_4jobs.sh
```

确认后提交：

```bash
CONFIRM_SUBMIT_RNA_EVIDENCE=1 \
bash 03_scripts/05_downstream_analysis/05_rna_umap/115_submit_rna_expression_evidence_4jobs.sh
```

4个作业全部完成后提交汇总：

```bash
csub < 03_scripts/05_downstream_analysis/05_rna_umap/117_submit_rna_expression_evidence_summary_8core.csub
```

主要汇总输出：

```text
08_qc_summary/35_rna_expression_evidence_audit/
  rna_expression_evidence_audit_summary.tsv
  rna_expression_evidence_per_group_rule.tsv
  rna_expression_evidence_per_well_rule.tsv
  rna_expression_evidence_per_well_validation.tsv

10_downstream_analysis/05_rna_expression_evidence_audit/
  rna_expression_evidence_audit.png
  rna_expression_evidence_audit.pdf
```

只有严格模式在GM12878-only孔中恢复human信号且仍保留足够细胞时，才将该模式合并成正式RNA表达矩阵并进行UMAP；否则RNA部分应报告为DNA残留影响严重、当前数据不足以可靠完成表达聚类。

## 步骤113–117结果：exon模式可作低深度挽救，不能作正式表达聚类

96个BAM及计数不变量均通过。三种模式的总exact UMI为：

- GeneFull：416,817；
- exon reverse：16,235；
- junction GeneFull：1,563。

在`UMI>=5且gene>=3`时，exon reverse保留334个细胞、5,766个UMI。其中前48个GM12878-only孔保留72个细胞，71个human、0个pig、1个tie；human exact-UMI加权比例为97.09%。说明exon证据显著改善物种方向，但细胞数和深度不足。

提高到`UMI>=10且gene>=5`后仅剩108个细胞；junction模式在最低阈值也仅剩12个细胞。因此：

- junction模式不进入矩阵；
- exon模式建立完整raw/filtered稀疏矩阵，保留全部16,235个原始exon UMI；
- 334个过滤后细胞仅用于探索性UMAP和给同一`物理孔|CB`的DNA细胞提供候选RNA物种标签；
- 不能表述为完整或高质量RNA表达聚类。

步骤118–121运行：

```bash
# 建立不覆盖旧GeneFull结果的完整exon salvage矩阵
csub < 03_scripts/05_downstream_analysis/05_rna_umap/119_submit_rna_exon_salvage_matrix_16core.csub

# 上一步validated后，运行探索性RNA UMAP
csub < 03_scripts/05_downstream_analysis/05_rna_umap/121_submit_rna_exon_salvage_umap_8core.csub
```

矩阵输出：

```text
10_downstream_analysis/06_rna_exon_salvage_matrix/
  raw_feature_bc_matrix/
  filtered_feature_bc_matrix_minUMI5_minGene3/
  rna_exon_salvage_complete_sparse.h5
  rna_exon_salvage_cell_metadata.tsv.gz

08_qc_summary/36_rna_exon_salvage_matrix/
  rna_exon_salvage_matrix_summary.tsv
  rna_exon_salvage_matrix_per_well.tsv
```

探索性UMAP输出：

```text
10_downstream_analysis/07_rna_exon_salvage_umap/
  rna_exon_salvage_umap_coordinates.tsv.gz
  rna_exon_salvage_selected_features.tsv.gz
  rna_exon_salvage_umap.png
  rna_exon_salvage_umap.pdf

08_qc_summary/37_rna_exon_salvage_umap/
  rna_exon_salvage_umap_summary.tsv
```

RNA UMAP按human/pig平衡选择重复检出的可变基因，但human和pig使用不同参考基因ID，因此物种分离仍部分来自参考特征空间，属于内部一致性展示，不是独立生物学验证。

## 步骤118–121结果：RNA探索性UMAP仍受极低深度限制

`minUMI5_minGene3`共保留334个细胞，来自86个物理孔；其中human 186、pig 146、tie 2。过滤后细胞的exon exact UMI中位数为13.5，基因数中位数为4；259/334个细胞检测到的基因数不超过5，126/334个细胞的exon exact UMI不超过10。矩阵为334 × 2,204，非零比例仅0.218%。

在原始稀疏表达空间中，10近邻同物种比例为67.1%，仅有限高于多数物种基线56.0%；UMAP空间中的10近邻同物种比例为70.4%。约73.7%的细胞对没有任何共同的非零基因。图中的细碎轨迹和小岛主要反映极稀疏离散计数，不能解释为稳定的转录细胞群。

UMAP实际选入45个特征（human 18、pig 27），包括多个物种特异参考ID及线粒体基因。因而有限的人猪邻域一致性同时包含参考特征空间效应，不能作为独立的GM12878/PK15验证。

图右下方约54个细胞形成的紧密“小岛”不是可解释细胞群：其中52个细胞在45个入选特征上完全没有非零值，整个334细胞中也恰好有52个此类细胞。该岛是零/近零特征向量被降维后堆叠形成的技术结构，应在图注中明确，不进行marker或群体解释。

RNA阶段结论：

- 完整exon salvage raw/filtered矩阵保留，作为可追溯输出；
- 当前UMAP仅作为低深度QC和方法探索图，不作为正式表达分群；
- 与DNA连接时，主分析标签使用`exon UMI>=10、gene>=5、major species fraction>=0.90`的66个高置信候选细胞；`UMI>=5、gene>=3、major fraction>=0.90`的200个细胞只用于敏感性分析；
- 后续DNA UMAP继续保持DNA特征独立拟合，RNA标签只在拟合完成后通过`物理孔|CB`连接、着色和计算一致率。
