# snMAR-seq DNA/RNA analysis pipeline

本仓库保存 2607 snMAR-seq 人/猪混合细胞文库中**实际采用并通过结果校验**的分析代码，用于后续追溯、修改和复现。原始测序数据、BAM、矩阵、BigWig、质控结果和作业日志不纳入 Git 管理。

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
- 人和猪的 contig 均带物种前缀，在 mapping、去重、矩阵和 BigWig 阶段保持分离。

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
    ├── 01_input_validation/     全量输入校验与任务清单
    ├── 02_dna_deduplication/    DNA CB-aware PCR 去重
    ├── 03_rna_expression/       RNA GeneFull UMI 矩阵
    ├── 04_dna_cpg_qc/           DNA CpG 汇总质控
    ├── 05_dna_single_cytosine/  DNA 单胞嘧啶调用与样本校验
    ├── 06_dna_signal_tracks/    WCG/GCH pseudo-bulk 轨迹与 BigWig
    ├── 07_matrix_packaging/     DNA/RNA 完整稀疏矩阵封装
    ├── 08_single_cell_qc/       完成度和单细胞覆盖质控
│   └── 90_optional_validation/  曾用于方法决策的可选 Probe
└── 05_downstream_analysis/      100-kb/TSS/ortholog 矩阵、QC 审计与 UMAP
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

当前整理版本结果：31 项通过；1 项 HDF5 测试在本机缺少 `h5py` 时跳过，服务器 `snTAPS_align` 环境中应执行。

## 版本范围

本仓库从本地 `snmar_2607_project_code_v19` 整理而来。旧版错误代码、一次性删除脚本及未采用的组合入口未迁入；原始 v19 备份没有被修改。

2026-08-31 的迁移分支新增了旧服务器最后导出的完整代码快照、尚未收尾的
`05_downstream_analysis` 阶段，以及物种分开的 IGV BigWig 入口。快照位于
`archive/e-drive-import-20260831/`，与整理后的主线并存，避免丢失正在进行的
版本。迁移时的数据保留边界见
[`docs/DATA_MIGRATION_20260831.md`](docs/DATA_MIGRATION_20260831.md)。
