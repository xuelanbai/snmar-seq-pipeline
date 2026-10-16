# Changelog

## v2026.10.16 — 2026-10-16

本版本是 2607 snMAR-seq 人/猪混合细胞项目的代码整理发布快照。

### 正式纳入

- DNA/RNA 严格 barcode 拆分、R1-only RNA fastp/STAR 和 paired-end DNA BWA 流程；
- CB-aware DNA PCR 去重与 reverse GeneFull exact-UMI RNA 矩阵；
- WCG、GCH、GCG 分离定义下的单胞嘧啶矩阵、单细胞 QC 和物种分离 BigWig；
- 100 kb 单细胞矩阵、96 孔汇总、跨孔 CB 审计、细胞阈值审计及 DNA 100 kb UMAP；
- IGV 单物种参考可读的去前缀 BigWig 输出。

### 探索性代码

- RNA exon-supported salvage 与低深度 UMAP；
- TSS±1.5 kb promoter UMAP、高变 promoter 敏感性分析；
- human/pig 同源 100 kb 特征和 Seurat 输入；
- reference-aware promoter 二维对照。

探索性代码随版本保留，但不标记为已获得最终生物学验证。

### 质量检查

- 本地 `unittest`：96 项通过，1 项因需要服务器 `h5py` 环境而跳过；
- Python 源文件通过编译检查；
- Git diff whitespace 检查通过；
- 未发现 GitHub token、AWS access key 或私钥文本。

### 兼容性与限制

- 集群提交脚本面向当前 `csub` 调度环境；
- 部分脚本保留 2607 项目路径作为兼容默认值，新项目应使用 `config/project.env` 覆盖；
- 100 kb reference-aware UMAP 不能单独证明仅由甲基化状态识别物种；
- RNA 数据受低拆分率、低深度与非特异性捕获影响，salvage UMAP 仅作探索。
