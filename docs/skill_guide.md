---
AIGC:
    Label: "1"
    ContentProducer: 001191440300708461136T1XGW3
    ProduceID: e82e471de5dc79d6e88eb5b8f99c8085_5f5ce68caeac11f1ac01525400e6dd8f
    ReservedCode1: nYv49zwuj6X8v46gInE3cu/8WNpSNu3j+nt53Md7xmrb7l0vFtnVi0MgN3VOOmpTR5QwlOcUo87EMwYxPuSwIWVQu0xtxRnxKM3zJ8YFiN769ZwbFze9JbwUSo1coyVkFkgtvIDxswi3yvzZ0JcxEH+uaZEq/+5vT+csKZ7lPqR4a7njNWw4VGY06JI=
    ContentPropagator: 001191440300708461136T1XGW3
    PropagateID: e82e471de5dc79d6e88eb5b8f99c8085_5f5ce68caeac11f1ac01525400e6dd8f
    ReservedCode2: nYv49zwuj6X8v46gInE3cu/8WNpSNu3j+nt53Md7xmrb7l0vFtnVi0MgN3VOOmpTR5QwlOcUo87EMwYxPuSwIWVQu0xtxRnxKM3zJ8YFiN769ZwbFze9JbwUSo1coyVkFkgtvIDxswi3yvzZ0JcxEH+uaZEq/+5vT+csKZ7lPqR4a7njNWw4VGY06JI=
---

# Skill: glycan-etl（糖类数据库 PDF→结构化入库 批量流水线）

> 状态：可用 | 版本 v1.0（2026-09-12）| 基于 glycan_etl_v3.py v3.1
> 适用操作系统：macOS（Windows 需将下述 python 改为对应 Python 路径，参数不变）

## 技能目标

把「糖类文献补充材料 PDF（Supporting Information，含 ¹H/¹³C NMR、2D 谱、结构式）」
自动解析为结构化记录并写入 PostgreSQL 糖类数据库（glycan_db）。
技能包提供**一键批量执行**：扫描目录下全部 PDF，逐个 dry-run 解析 → 自动生成汇总报告 →
可选正式入库。避免每次手动改配置、逐个跑命令。

## 技能位置与文件清单

仓库根目录即技能包根目录（`glycan-db/`）。

| 文件 | 说明 |
|------|------|
| `glycan_etl/core.py` | 核心解析入库引擎（正文行内式：单糖/寡糖/多糖 + 1D/2D 谱） |
| `glycan_etl/table_parser.py` | 表格式文献解析器（分子量表 / 位移归属表 / 正文键连） |
| `glycan_etl/embeddings.py` | BGE-M3 向量化（回写 `nmr_shifts_1d.embedding`） |
| `pipelines/batch_etl.py` | 批量执行器（本技能主入口，扫描目录/多文件批量处理） |
| `pipelines/config.yaml` | 配置示例（文献溯源 + 数据库连接），真实配置另存 `config.local.yaml` |
| `etl_reports/` | 批量执行自动生成的报告目录（CSV + Markdown，时间戳命名） |
| `db/init.sql` | 一键幂等建库（含向量列）；分步版见 `db/migrations/v001…v004` |
| `db/schema_design.md` | schema 设计文档 |

## 数据库表

`sugars` / `residues` / `nmr_experiments` / `nmr_shifts_1d` / `nmr_correlations_2d` /
`physicochemical` / `polysaccharide_props` / `literature`
（建表二选一：① 一键 `db/init.sql`；② 分步 `db/migrations/v001_mvp.sql` →
`v002_extend.sql` → `v003_qc.sql` → `v004_vector_and_provenance.sql`。
向量检索列 `nmr_shifts_1d.embedding` 由 v004 创建，漏掉它 `embeddings.py` 会报列不存在）

## 触发场景

- 用户说：把这份/这些 SI.pdf 入库、跑 ETL、解析糖类 NMR 数据、批量处理补充材料
- 用户提供一份或多份 PDF（或指一个包含 PDF 的目录）

## 执行步骤（标准流程）

### 0. 环境检查
```bash
python3 -c "import pdfplumber, psycopg2, yaml"   # 缺哪个: pip install 对应包
```

### 1. 自检（确认环境与引擎正常）
```bash
python3 pipelines/batch_etl.py --self-test
```
输出 `SELF-TEST PASS` 即正常。

### 2. 批量 dry-run（必做，零风险）
```bash
# 处理一个目录下所有 PDF（可用 --filter 只挑含某关键词的文件）
python3 pipelines/batch_etl.py --dir "/path/to/补充材料目录" [--filter SI]

# 或显式指定多份 PDF
python3 pipelines/batch_etl.py --pdf a.pdf b.pdf c.pdf
```
生成报告 `etl_reports/manifest_<时间戳>.md`（含：每份识别记录数 / 类型分布 / QC 状态 / 明细）。
识别为 0 条的 PDF：先查文本层质量或结构写法，再决定是否入库。

### 3. 准备入库配置
编辑 `etl_config.yaml`：
- `db.*`：改成真实 PostgreSQL 连接（host/port/dbname/user/password）
- `doi / journal / year`：批量多份文献时更推荐用 `meta.yaml` 逐文件指定（见下）

示例 `meta.yaml`（可选，按文件名逐份指定文献溯源）：
```yaml
fmicb-11-00806.pdf:
  doi: 10.3389/fmicb.2020.00806
  journal: Front. Microbiol.
  year: 2020
dkw377.pdf:
  doi: 10.1093/.../dkw377
  journal: ...
  year: ...
```

### 4. 批量正式入库
```bash
# 单套文献元数据（同一篇文献的多分 SI 等）:
python3 pipelines/batch_etl.py --dir "/path/to/目录" --config pipelines/config.local.yaml

# 多份不同文献（推荐，用 meta 映射）:
python3 pipelines/batch_etl.py --dir "/path/to/目录" --config pipelines/config.local.yaml --meta meta.yaml
```
日志末尾出现「入库 N 条」即成功。**QC flagged 的记录默认仍入库并带 `qc_status='flagged'`
标记**（便于回溯）；确需丢弃时加 `--skip-flagged`。

### 5. 给用户汇报
- 报告路径（Markdown + CSV）
- 每份 PDF 识别条数、flagged 条数与原因、入库条数
- 识别为 0 条的 PDF 单独点名，说明可能原因与建议

## 注意事项 / 避坑

1. **先 dry-run 再入库**：严禁跳过 dry-run 直接连库，先查识别率与溶剂/温度正确性。
2. **QC flagged 默认保留**：质控标记（R1 位移越界 / R2 J 耦合异常 / R3 构型-J 不符等）
   的记录会入库并保留 `qc_status='flagged'`，下游按状态过滤即可；确需丢弃加 `--skip-flagged`。
3. **重复导入安全**：同一 PDF 重跑不会重复累积（按 结构+文献+谱类型+溶剂 先清后写）。
4. **结构编码缺失不会串档**：解析不出 GlycoCT 时用「结构指纹+DOI」生成占位编码，
   不同结构彼此隔离，不会互相覆盖。
5. **溯源元数据**：批量多文献务必用 `--meta`，否则全部 PDF 会共用 config 的同一 doi/journal/year，导致溯源串档。
6. **解析范围**：仅支持正文含 `.pdf` 文本层的 SI（扫描件需先 OCR）；结构写法仅支持 IUPAC 连接式 /
   单糖标题行 / 多糖重复单元写法（`→4)-β-D-Glcp-(1→` 与 `(1→4)-linked π-D-Glcp`）。
7. **依赖**：`pip install pdfplumber psycopg2-binary pyyaml`
8. 中间 JSON（完整解析明细）可用 `--save-json` 开启，与报告同目录输出。

## 常见问题

| 现象 | 处理 |
|------|------|
| `SELF-TEST` 失败 | 重新 `pip install pdfplumber psycopg2-binary pyyaml` 后重试 |
| 连不上库 | 检查 `pipelines/config.local.yaml` 的 host/port/dbname/user/password |
| 某份 PDF 0 条记录 | 无文本层(需OCR) 或结构写法不支持，看 dry-run 报告定位 |
| `--config` 文件缺失 | batch 脚本自动降级为 dry-run，不会误写库 |
*（内容由AI生成，仅供参考）*
