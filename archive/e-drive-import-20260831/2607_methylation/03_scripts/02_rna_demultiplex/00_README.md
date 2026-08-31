# 02 RNA 严格10 bp I7拆分

本目录只处理本批snMAR RNA文库：

```text
R1 = transcript read
R2 = T7(8) + RNA-linker(22–25) + T5(8) + UMI(8) + dT(31) + read
FASTQ header index = I7 + I5
cell barcode = T7 + T5 + I7 + I5
```

## 正式主分析口径

- R2前8 bp必须精确匹配T7；
- 随后的linker必须精确匹配该T7的RNA-linker；
- linker后8 bp必须精确匹配T5；
- I5必须完整10 bp精确匹配；
- I7必须是实际观测到的完整10 bp，并精确匹配白名单；
- UMI只提取、不纠错；dT完整度只做QC；
- `CTCCGACGC` 后插入 `GTCA` 的linker属于DNA文库残留，只计数并丢弃；
- 8 bp I7即使对应唯一白名单前缀，也缺少两个实际观测碱基，只作为被拒绝的
  潜在救援候选计数，不补全、不写FASTQ。

正式输出全部是完整10 bp I7：

```text
CR:Z:原始T7+T5+10bp-I7+I5
CB:Z:严格匹配后的T7+T5+10bp-I7+I5
CN:Z:T7_T5_I7_I5名称
UR:Z:原始UMI
UB:Z:UMI（当前不纠错）
XD:i:0
XE:Z:0,0,0,0
XL:i:RNA-linker长度
XI:Z:full10
```

## 为什么需要重新统计

旧 `full_strict` 统计把唯一8 bp I7前缀补全后纳入，得到42.670024%的分配率。
该结果只保留为预分析记录，不是最终严格结果。

新代码使用独立运行名，避免误用旧分片：

```text
full_strict10
write_strict10
```

新统计同时报告：

```text
I7_full10_exact
I7_prefix8_unique_candidate_rejected
I7_rejected
I7_unmatched
potential_prefix8_rescue_structure_I5_exact
strict_assigned_I7_full10
```

其中 `potential_prefix8_rescue_structure_I5_exact` 只量化“如果未来有充分证据，理论上
可救回多少”，本轮仍全部拒绝。

## 执行顺序

```bash
cd /share/org/YZWL/yzwl_wujx/Projects/2607_methylation

bash 03_scripts/02_rna_demultiplex/02_setup_project_dirs.sh
python3 03_scripts/02_rna_demultiplex/07_verify_install.py

# 重新跑只接受完整10 bp I7的全量统计
csub < 03_scripts/02_rna_demultiplex/09_submit_full_stats_64core.csub
```

检查：

```text
08_qc_summary/02_rna_demultiplex/full_strict10_summary.tsv
08_qc_summary/02_rna_demultiplex/full_strict10_per_pair.tsv
08_qc_summary/02_rna_demultiplex/full_strict10_top_unknown_barcodes.tsv
```

确认最终严格分配量足以支持下游分析后，再提交：

```bash
csub < 03_scripts/02_rna_demultiplex/10_submit_split_exact_32core.csub
```

正式输出：

```text
04_fastq_split/rna_strict10/<sample>/*.assigned_R1.fq.gz
04_fastq_split/rna_strict10/<sample>/*.assigned_R2.fq.gz
08_qc_summary/02_rna_demultiplex/write_strict10_summary.tsv
08_qc_summary/02_rna_demultiplex/write_strict10_per_pair.tsv
```

步骤11/12保留为已完成的只读结构诊断。步骤13/14因GTCA插入型已确认为DNA文库
残留而取消。
