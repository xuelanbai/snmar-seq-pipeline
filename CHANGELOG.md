# Changelog

## 1.0.0 - 2026-08-03

- 以 `snmar_2607_project_code_v19` 为基线建立可追溯仓库。
- 按功能细分 DNA/RNA fastp、reference、mapping 和 post-mapping 目录。
- 保留最终使用的 RNA R1-only STAR、DNA CB-aware dedup、reverse GeneFull 和 WCG/GCH 单碱基流程。
- 将方法选择 Probe 隔离到 optional validation 目录。
- 排除旧版错误逻辑、一次性清理脚本和未采用的组合入口。
- 新增 `PROJECT`、`SNMAR_SCRIPT_ROOT`、`CONDA_SH` 等可覆盖配置。
- 新增中文逐脚本输入输出说明、版本删减记录、数据忽略规则和测试说明。
- 本地单元测试：31 passed，1 skipped（服务器 h5py 环境执行）。
