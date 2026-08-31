# 01 DNA 拆分：严格率、救援率与正式输出

## 结构与 reference

```text
R1 前 8 bp = T5
R2 前 8 bp = T7
FASTQ header index = I7 + I5
cell barcode = T7 + T5 + I7 + I5
```

本批 T7-1..T7-96 与旧参考表前 96 条完全一致。旧参考脚本通过完整 barcode
字典查找，默认只支持完全匹配；没有汉明距离或错配参数。

## 当前正式结果

已完成 768 对 FASTQ 的 `write_rescue`：

```text
total_pairs       19,970,035,631
strict/XD0        17,221,376,297（86.236082%）
rescued/XD1        1,674,949,100（ 8.387312%）
rescued/XD2          140,924,698（ 0.705681%）
assigned          19,037,250,095（95.329074%）
unassigned           932,785,536（ 4.670926%）
```

正式文件位于：

```text
04_fastq_split/write_rescue/
08_qc_summary/01_demultiplex/write_rescue_*
09_logs/01_demultiplex/write_rescue/
```

这些文件保留。后续以 `XD=0` 为主分析，`XD=1/2` 只做敏感性分析；不再提交
`11_submit_split_exact_32core.csub`，避免重读约 200 亿对 reads 并重复占用磁盘。

## 实际拆分算法

1. `03_discover_fastq_pairs.py` 递归发现 R1/R2，排除 `Undetermined`；
2. `04_barcode_qc_demultiplex.py` 同步读取配对 FASTQ；
3. 提取 R1-T5、R2-T7、header-I7/I5；
4. 四段完全匹配记为 `XD=0`；
5. `write_rescue` 中，每段仅允许唯一 1-bp 最近邻，四段总错配不超过 2；
6. 歧义候选、超过阈值、短 read 或 header 无法解析均不分配；
7. 原序列不修剪，输出仍按原始 pair 分片；
8. header 写入：

```text
CR:Z:原始 T7+T5+I7+I5
CB:Z:匹配后的 barcode 序列
CN:Z:T7_T5_I7_I5 名称
XD:i:总错配数
XE:Z:T7,T5,I7,I5 各段错配数
```

v5 若重新写 FASTQ，会把 comment 规范为 SAM 可继承的 tab 分隔标签。服务器现有
`write_rescue` 是上一版已完成结果，不需要重跑；Mapping 时用
`../03_mapping/01_stream_strict_fastq.py` 流式筛选 `XD=0`，同时规范旧 header。

## 已完成与可选命令

以下 Probe、全量统计和正式 rescue 已完成，不应重复提交：

```text
09_submit_probe_32core.csub
10_submit_full_stats_64core.csub
12_submit_split_rescue_32core.csub
```

`13_cleanup_test_intermediates.sh` 默认 dry-run，只清理已汇总的 Probe/full-stats
分片和日志，不触碰原始数据、汇总表或 `write_rescue`。

```bash
# 先查看
bash 03_scripts/01_demultiplex/13_cleanup_test_intermediates.sh

# 只有确认六张 Probe/full-stats 汇总 TSV 已备份后才执行
CONFIRM_CLEANUP=YES \
bash 03_scripts/01_demultiplex/13_cleanup_test_intermediates.sh
```
