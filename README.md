# snMAR-seq DNA/RNA analysis pipeline

> **代码整理发布版：`v2026.10.16`**<br>
> 对应 2607 snMAR-seq 人/猪混合细胞文库；发布日期：2026-10-16。

本仓库保存本批数据从 FASTQ 拆分、mapping、去重、单碱基 WCG/GCH 矩阵到下游降维分析的可追溯代码。原始测序数据、BAM、矩阵、BigWig、质控结果和作业日志不纳入 Git 管理。生产主流程、已运行的探索性分析和仅代码就绪的方法试验在下表中分别标注，不能仅因脚本存在就将其视为已获得生物学验证。

## 发布状态

| 范围 | 状态 | 说明 |
|---|---|---|
| 步骤 01–86 | **正式主流程（已完成）** | 严格拆分、fastp、联合参考 mapping、CB-aware DNA 去重、RNA exact-UMI、单碱基矩阵、QC 与 BigWig |
| 步骤 87–110 | **正式下游（已完成）** | 100 kb 矩阵、96 孔汇总/CB 审计、细胞 QC、DNA 100 kb UMAP 与 GM12878/PK15 内部一致性审计 |
| 步骤 111–121 | **探索性（已运行）** | RNA exon-supported salvage；受文库深度和非特异性捕获影响，不作为独立可靠的细胞类型结论 |
| 步骤 122–134 | **探索性（已运行）** | TSS±1.5 kb promoter UMAP；本批结果未形成清晰分群，保留用于方法追溯 |
| 步骤 135–146 | **方法试验（非最终结果）** | human/pig 同源 100 kb 特征与 Seurat 输入流程；不得与已验证主流程混称 |
| 步骤 147–154 | **敏感性分析（非最终结果）** | 高变 promoter 与 reference-aware promoter 对照；用于定位稀疏度和参考序列差异影响 |

版本、验证范围和已知限制见 [v2026.10.16 发布说明](docs/发布说明_v2026.10.16.md)；逐文件输入输出见 [分析流程与代码说明](docs/分析流程与代码说明.md)。

## 最终采用的分析方案

- DNA 主分析仅使用严格条形码匹配 reads（`XD:i:0`）。
- RNA 使用完整 10 bp I7 严格匹配结果（`XD:i:0, XI:Z:full10`）。
- RNA 在提取 CB/UMI 后仅保留 R1；R2 不进入 fastp 和 mapping。
- RNA 采用单端 STAR，比对到带物种前缀的人+猪联合参考基因组。
- DNA 采用 paired-end BWA-MEM，比对到同一联合参考基因组。
- DNA PCR 去重在 CB 内进行，键为 `CB + reference + unclipped template endpoints + orientation`。
- RNA 采用 reverse GeneFull，按 `(sample, CB, gene, exact UMI)` 去重计数。
- DNA 单碱基矩阵区分：
  - `WCG`：内源甲基化信号；
  - `GCH`：GpC 酶标记的开放性信号；
  - `GCG`：CpG/GpC 重叠歧义位点，仅审计，不进入两类主信号。
- 人和猪的 contig 均带物种前缀，在 mapping、去重、矩阵和主 BigWig 阶段保持分离；另由 step 85 输出只移除物种前缀、供单物种参考在 IGV 中展示的 BigWig。

实验单位说明：本批数据包含 **96 个物理孔**，后48孔加入猪细胞，用于检验物种分离和 UMAP；这些孔是否来自同一份生物学样本需由实验记录确认，流程不自行推断。上游为兼容已有文件名保留的 `sample` 字段实际表示物理孔。跨孔审计确认本批36 bp CB没有任何重复，因此原统计值不变；内部仍保留 `well_index|CB` 便于来源追踪。总体信号和 BigWig 可跨96孔汇总，但必须按 human/pig 参考物种分开。DNA UMAP 不按物种预拆分，物种标签仅在 UMAP 完成后用于着色和准确性验证；DNA/RNA 通过 CB 匹配并核对孔号。

## 目录

```text
03_scripts/
├── 01_demultiplex/              DNA 条形码质控和严格拆分
├── 02_rna_demultiplex/          RNA strict10 条形码质控和拆分
├── 03_mapping/
│   ├── 01_dna_fastp/            DNA fastp
│   ├── 02_rna_fastp/            RNA R1-only fastp
│   ├── 03_joint_reference/      人+猪联合参考及索引
│   ├── 04_rna_star/             RNA 单端 STAR
│   └── 05_dna_bwa/              DNA paired-end BWA
├── 04_postmapping/
│   ├── 01_input_validation/     全量输入校验与任务清单
│   ├── 02_dna_deduplication/    DNA CB-aware PCR 去重
│   ├── 03_rna_expression/       RNA GeneFull UMI 矩阵
│   ├── 04_dna_cpg_qc/           DNA CpG 汇总质控
│   ├── 05_dna_single_cytosine/  DNA 单胞嘧啶调用与样本校验
│   ├── 06_dna_signal_tracks/    WCG/GCH pseudo-bulk 轨迹与 BigWig
│   ├── 07_matrix_packaging/     DNA/RNA 完整稀疏矩阵封装
│   ├── 08_single_cell_qc/       完成度和单细胞覆盖质控
│   └── 90_optional_validation/  曾用于方法决策的可选 Probe
└── 05_downstream_analysis/
    ├── 01_dna_100kb_bins/       DNA 单细胞 100 kb WCG/GCH 稀疏矩阵
    ├── 02_cross_well_summary/   96孔、分物种的 DNA/RNA 汇总与 CB 审计
    ├── 03_cell_qc_filtering/    细胞阈值、跨模态一致性及RNA物种/残留审计
    ├── 04_dna_umap/             DNA 100 kb WCG、GCH和联合UMAP
    ├── 05_rna_umap/             RNA exon-supported探索性矩阵与UMAP
    ├── 06_dna_tss_umap/         promoter/TSS UMAP及敏感性分析
    └── 07_dna_ortholog_100kb/   human/pig同源100 kb方法试验
```

完整脚本索引、逐步输入和输出见 [分析流程与代码说明](docs/分析流程与代码说明.md)。删除和保留依据见 [版本整理与删减记录](docs/版本整理与删减记录.md)。

## 配置

```bash
git clone https://github.com/xuelanbai/snmar-seq-pipeline.git
cd snmar-seq-pipeline
cp config/project.env.example config/project.env
vim config/project.env
source config/project.env
```

`PROJECT` 指向数据与结果根目录，`SNMAR_SCRIPT_ROOT` 指向本仓库的 `03_scripts`。这样仓库不必克隆到数据目录中。

当前脚本仍保留本项目服务器路径作为兼容默认值，以保证对已经完成的 2607 项目可以精确追溯；新项目应通过 `config/project.env` 或命令行参数覆盖路径。

## 运行原则

1. `.csub` 文件用于当前集群的 `csub` 调度器。
2. 名称包含 `template` 的 `.csub` 不能直接提交，必须由相邻 launcher 渲染任务编号。
3. launcher 默认 dry-run；只有设置对应的 `CONFIRM_SUBMIT_*=1` 才会提交。
4. 每个全量阶段完成后先运行 summary/validation，再进入下一阶段。
5. 单碱基矩阵中的缺失表示无覆盖（`NA`），不能解释为 0% 信号。

## 本地检查

无需真实 FASTQ/BAM 的单元测试：

```bash
python -m unittest discover -s tests -v
```

`v2026.10.16` 在本地通过 96 项测试；另有 1 项依赖 `h5py` 的 HDF5 测试按设计在服务器 `snTAPS_align` 环境验证。测试通过表示代码级计数守恒、边界条件和提交入口符合预期，不等同于所有探索性分析均获得生物学阳性结论。

## RNA 物种与DNA残留审计（步骤 98–102）

RNA严格拆分统计已经单独记录 `DNA_residual_GTCA_linker_given_T7`；该指标描述的是原始RNA FASTQ中被识别为DNA结构、因而没有进入严格RNA输出的reads，不能单独证明DNA残留已经泄漏到RNA矩阵。步骤98–102在不修改BAM和矩阵的前提下进一步检查：

- 4个并行作业扫描96个RNA BAM，区分primary、`NH:i:1`唯一比对、human/pig contig以及剪接/未剪接reads；
- 汇总时连接拆分统计、GeneFull JSON、exact UMI和两个候选细胞阈值；
- 前48孔按实验设计标为GM12878，仅用于检验RNA物种结果，不作为RNA物种判定本身；
- `validated`只表示统计守恒，最终的 `biological_species_validation` 才表示纯人孔是否符合预期。

提交前先预览：

```bash
bash 03_scripts/05_downstream_analysis/03_cell_qc_filtering/100_submit_rna_species_alignment_audit_4jobs.sh
```

确认后提交4个32核作业：

```bash
CONFIRM_SUBMIT_RNA_SPECIES_AUDIT=1 \
bash 03_scripts/05_downstream_analysis/03_cell_qc_filtering/100_submit_rna_species_alignment_audit_4jobs.sh
```

4个作业全部完成后提交步骤102汇总。

## DNA 100 kb UMAP（步骤 103–108）

主分析使用 `D_total_20k`（4,778个细胞），并行运行
`D_total_50k`（2,994个细胞）作为严格敏感性分析。WCG、GCH和两者联合
分别输出UMAP。物种和孔组不参与细胞过滤、特征选择、SVD或UMAP，只在降维
完成后着色并统计近邻一致率。无覆盖bin按NA处理为bin总体期望的零残差，
不解释为0%甲基化或0%开放性。

先预览两个32核全量作业：

```bash
bash 03_scripts/05_downstream_analysis/04_dna_umap/106_submit_dna_100kb_umap_2jobs.sh
```

确认后并行提交：

```bash
CONFIRM_SUBMIT_DNA_UMAP=1 \
bash 03_scripts/05_downstream_analysis/04_dna_umap/106_submit_dna_100kb_umap_2jobs.sh
```

两项均完成后提交
`108_submit_dna_100kb_umap_summary_4core.csub`，比较D20与D50的稳定性。

## 版本范围

本仓库从本地 `snmar_2607_project_code_v19` 整理而来。旧版错误代码、一次性删除脚本及未采用的组合入口未迁入；原始 v19 备份没有被修改。

## DNA TSS±1.5 kb UMAP（步骤 122–134）

100 kb UMAP 后新增了共享 promoter 特征分析：先建立经染色体边界截断、
互不重叠的 TSS±1.5 kb promoter cluster，再将 human/pig 的保守同名
单基因 promoter 投影到共享轴，分别拟合 WCG、GCH 和 joint UMAP。
GM12878/PK15 mapping 标签仅在拟合完成后用于着色与内部一致性检查。
完整提交顺序、输入输出和解释边界见
[DNA TSS UMAP 步骤 122–134](docs/DNA_TSS_UMAP_步骤122-134.md)。

## 后续探索性方法

- RNA exon-supported salvage（111–121）：见主追溯文档；输出仅用于判断低深度数据是否仍保留表达结构。
- human/pig 同源 100 kb（135–146）：见 [DNA同源100kb Seurat试验](docs/DNA同源100kb_Seurat试验_步骤135-146.md)。
- 高变 promoter（147–150）：见 [DNA TSS高变promoter](docs/DNA_TSS高变promoter_步骤147-150.md)。
- 非同源 promoter / 同源 100 kb 对照（151–154）：见 [DNA二维对照](docs/DNA二维对照_非同源promoter与同源100kb_步骤151-154.md)。

这些方法保留在发布版中用于复现方法尝试和继续开发，但不改变本批数据目前的核心结论：100 kb 的物种分离包含参考序列差异贡献，promoter 与 RNA 探索性结果受覆盖稀疏和文库质量限制。
