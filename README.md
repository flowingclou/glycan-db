---
AIGC:
    Label: "1"
    ContentProducer: 001191440300708461136T1XGW3
    ProduceID: e82e471de5dc79d6e88eb5b8f99c8085_5d581970af3411f1ac01525400e6dd8f
    ReservedCode1: 85y0F5zhIRxYR8A9oOV3a7DKtpZUkFFOUfD+78b+dgyW5bxhFe+5TzdLOSlZQeR7Q2CN/cWVozO95z6QJUncmbhA7vAbPXXJFJsf8gkvRika90U9D5RWw4DTM2yi8NQztvvKgoEMxNkHmon+RqzHe4F0SlqF7nWZV5MqqNZRl2mlPMr4sjaQzfY4C6Y=
    ContentPropagator: 001191440300708461136T1XGW3
    PropagateID: e82e471de5dc79d6e88eb5b8f99c8085_5d581970af3411f1ac01525400e6dd8f
    ReservedCode2: 85y0F5zhIRxYR8A9oOV3a7DKtpZUkFFOUfD+78b+dgyW5bxhFe+5TzdLOSlZQeR7Q2CN/cWVozO95z6QJUncmbhA7vAbPXXJFJsf8gkvRika90U9D5RWw4DTM2yi8NQztvvKgoEMxNkHmon+RqzHe4F0SlqF7nWZV5MqqNZRl2mlPMr4sjaQzfY4C6Y=
---

# glycan-db · 糖类结构数据库数据工程流水线

> 为 **AI 辅助推断多糖结构平台** 提供「已知结构 — 谱图」真值参考数据库的数据工程流水线。
> 从糖类文献补充材料（Supporting Information，含 ¹H/¹³C NMR、2D 谱、结构式）的 PDF
> 自动解析为结构化记录，写入 PostgreSQL**扩展了 pgvector 向量检索**的糖类数据库。

本仓库是模块化、可分享的工程化整理版本：把散装的 ETL 脚本重组为**解析引擎包 + 批量管道 + 数据库迁移 + 一键环境**，
便于版本管理、协作开发，也便于后续对接检索知识库与 AI 推理平台。

---

## 1. 项目定位

| 层面 | 说明 |
|---|---|
| 数据来源 | 多糖/寡糖文献补充材料 PDF（正文文本层，含 NMR 归属表、分子量表、结构写法） |
| 解析产出 | sugars / residues / nmr_experiments / nmr_shifts_1d / nmr_correlations_2d / physicochemical / polysaccharide_props / literature |
| 存储 | PostgreSQL 16 + pgvector（embedding 列） |
| 下游用途 | 作为 AI 辅助推断多糖结构的**真值/参照数据库**：检索知识库召回相似结构，再交给 AI 平台做结构推断验证 |

---

## 2. 快速开始

### 2.1 一键起库（PostgreSQL 16 + pgvector）

```bash
docker compose up -d
# 首次启动会自动执行 db/init.sql（幂等建库 + 种子数据，含 vector/pgcrypto 扩展）
```

默认连接：`localhost:5432`，库 `glycan_db`，用户 `glycan`，密码 `glycan_dev`（**仅限本地开发，生产请修改**）。

### 2.2 安装 Python 依赖

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

### 2.3 自检

```bash
python3 tests/run_tests.py     # py_compile 全仓校验 + 回归测试 + batch_etl --self-test
# 或单跑
python3 pipelines/batch_etl.py --self-test   # 输出 SELF-TEST PASS 即正常
python3 tests/test_regressions.py            # 11 项防回退回归测试（不依赖数据库）
```

### 2.4 批量解析 PDF

```bash
# ① 批量 dry-run（零风险，推荐先跑）
python3 pipelines/batch_etl.py --dir "/path/to/补充材料目录" [--filter SI]

# ② 准备真实连接配置（复制示例，勿把真实口令提交到仓库）
cp pipelines/config.yaml pipelines/config.local.yaml   # 编辑 db.* 为真实值

# ③ 批量正式入库
python3 pipelines/batch_etl.py --dir "/path/to/目录" --config pipelines/config.local.yaml
```

多份不同文献批量入库时，推荐用 `--meta meta.yaml` 按文件名逐份指定 doi/journal/year，避免溯源串档。

#### 入库行为说明（重要）

- **重复导入是安全的**：同一份 PDF 反复入库不会重复累积数据。脚本按
  `(结构, 文献, 谱类型, 溶剂)` 先清理旧的派生行再重写，可放心重跑。
- **质控命中默认仍入库**：QC 判定为 `flagged` 的记录会写入库并保留 `qc_status='flagged'`
  标记（便于回溯与人工复核）。若确实要丢弃，加 `--skip-flagged`。
- **结构编码缺失不会串档**：文献未给出可用 GlycoCT 时，脚本用
  结构指纹 + DOI 生成确定性占位编码（`UNRESOLVED:...`），不同结构彼此隔离。
  代价是同一结构被多篇文献报道时会各占一条记录，待补全真实 GlycoCT 后再归并——
  这比"不同结构被错并成一条、谱图混在一起"安全得多。

---

## 3. 目录结构

```text
glycan-db/
├── README.md               # 本文件
├── LICENSE                 # MIT
├── pyproject.toml(可选)    # （预留）包元数据
├── requirements.txt        # Python 依赖：pdfplumber / psycopg2-binary / pyyaml
├── docker-compose.yml      # PostgreSQL16 + pgvector 一键环境（挂载 init.sql）
├── .gitignore              # 排除缓存 / 虚拟环境 / 真实密钥配置文件
├── db/
│   ├── init.sql            # 四份迁移合并后的幂等一键建库脚本（可重复执行）
│   ├── migrations/         # 分版本迁移 SQL
│   │   ├── v001_mvp.sql    #   单糖四表（sugars/literature/nmr_experiments/nmr_shifts_1d）
│   │   ├── v002_extend.sql #   残基/2D 相关峰/物化/多糖性质扩展
│   │   ├── v003_qc.sql     #   质控触发器 R1-R10 + 全库校验函数
│   │   ├── v004_vector_and_provenance.sql
│   │   │                   #   pgvector embedding 列 + HNSW 索引 + 溯源列
│   │   └── v005_structure_level.sql
│   │                       #   结构表达等级 + 残基组成式
│   └── schema_design.md    # 数据库 schema 设计文档
├── glycan_etl/             # 解析引擎包
│   ├── __init__.py         #   包初始化 + __version__
│   ├── core.py             #   单份 PDF 解析入库引擎（v3.x 全功能）
│   ├── table_parser.py     #   表格式文献解析器（分子量表/归属表/正文键连）
│   ├── glycoct.py          #   标准 GlycoCT 结构编码生成 + 校验
│   └── embeddings.py       #   BGE-M3 向量化脚本（API Key 仅从环境变量读取）
├── pipelines/
│   ├── batch_etl.py        # 批量执行器（自检/批量 dry-run/批量入库/汇总报告）
│   ├── backfill_glycoct.py # 存量记录结构编码回填/修正
│   └── config.yaml         # 配置示例（含数据库连接占位符与溯源字段）
├── tests/
│   ├── run_tests.py        # 自检入口（py_compile + 回归测试 + self-test）
│   └── test_regressions.py # 防回退回归测试（结构覆盖/命名/标题抽取）
└── docs/
    └── skill_guide.md      # glycan-etl 技能使用说明（可用即用指南）
```

---

## 4. 模块说明

### 4.1 解析引擎包 `glycan_etl/`

| 模块 | 职责 |
|---|---|
| `core.py` | 核心解析入库引擎：PDF 文本层抽取、单糖/残基/键连结构写法识别、NMR 位移提取、QC 判定、入库。复用数据类 `GlycanRecord / Residue / NMRExperiment / Peak1D / PolysaccharideProps` 与 `MONOSACCHARIDES` 常量 |
| `table_parser.py` | 补齐表格式文献的解析：分子量表（RT/Mp/Mw/Mn/单糖 mol%）、位移归属表（位置数字表头 与 `H1/C1` 成对表头两种版式）、正文甲基化/糖苷键连；并从正文结论句提取**主链连接式**还原完整序列。含双栏页面线性化与 `x_tolerance` 词粘连处理 |
| `embeddings.py` | 读取库内结构化记录，调用 SiliconFlow BGE-M3 API 生成向量并回写 embedding 列。API Key 仅从环境变量 `SILICONFLOW_API_KEY` 读取，代码内不含任何明文密钥 |

### 4.2 批量管道 `pipelines/`

`batch_etl.py` 封装「单 PDF → 手动改 config → dry-run → 入库」多步流程为一次命令：
- `--self-test`：环境自检（内置样例文本，不依赖 PDF 与数据库）
- `--dir / --pdf`：整目录或指定多份 PDF
- `--dry-run`（缺省）：零风险解析，生成汇总报告到 `etl_reports/`
- `--config`：正式入库连接配置；缺省自动降级为 dry-run
- `--meta`：按文件名逐份指定文献溯源
- `--save-json`：额外输出完整解析明细
- `--skip-flagged`：质控违规记录不入库（默认入库并标记，见 §2.4）

### 4.3 数据库 `db/`

- `init.sql` 为四份迁移（v001 MVP → v002 扩展 → v003 质控 → v004 向量与溯源）合并后的
  **幂等脚本**：`CREATE EXTENSION IF NOT EXISTS vector/pgcrypto` + 全部 `CREATE TABLE IF NOT EXISTS`
  + 种子数据 `ON CONFLICT / NOT EXISTS` 防重，**可重复执行不报错、不产生脏数据**（已实测连续执行 3 次）。
- 手工分步建库可依次执行：
  `migrations/v001_mvp.sql → v002_extend.sql → v003_qc.sql → v004_vector_and_provenance.sql`。
  四份迁移同样逐条带 `NOT EXISTS` 防重，可重复执行。
- **向量检索能力来自 v004**：`nmr_shifts_1d.embedding vector(1024)` 列与 HNSW 余弦索引
  定义在 v004 中。缺了这一步，`glycan_etl/embeddings.py` 会直接报
  `column "embedding" does not exist`，README §5 的近邻检索 SQL 也无法执行。

### 4.4 结构编码与表达等级（`glycan_etl/glycoct.py`）

`sugars.glycoct` 存放**标准 GlycoCT**（`RES`/`LIN` 分段 + 规范 basetype），
可被 glypy、GlyTouCan 等外部工具解析：

```text
RES
1b:a-dglc-HEX-1:5
2b:x-dglc-HEX-1:5
LIN
1:2o(4+1)1d          # 残基2 的 O4 ← 残基1 的异头碳
```

编码规则要点：

- 连接写作 `<受体残基>o(<受体位点>+1)<供体残基>d`，供体提供异头碳、受体提供羟基氧；
- 糖醛酸用 `|6:a`（如 GalA）、6-脱氧糖用 `|6:d`（如 Rha/Fuc）；
- N-乙酰氨基糖按规范以独立取代基残基 `2s:n-acetyl` 表达；
- 生成结果用 glypy 校验（未安装 glypy 时跳过校验，不影响生成）。

**不是所有记录都能给出完整结构。** GlycoCT 只能表达确定结构，而文献里
大量多糖只给出甲基化/组成信息、残基间连接顺序未知。硬生成编码等于伪造
数据，因此 `sugars.structure_level` 显式记录表达程度：

| structure_level | 含义 | glycoct 字段 |
|---|---|---|
| `complete` | 单糖/寡糖，连接明确 | 完整结构编码 |
| `repeat_unit` | 单一重复单元多糖 | 一个重复单元的编码 |
| `composition_only` | 仅残基组成，连接顺序未知 | 空；结构信息在 `composition`（如 `GalA7,Ara3,Gal3,Rha2,GlcA`） |

下游按此字段决定能否做结构级比对：`complete` / `repeat_unit` 可直接比对，
`composition_only` 只能按组成筛选。存量记录可用回填脚本修正：

```bash
python3 pipelines/backfill_glycoct.py --config pipelines/config.local.yaml --dry-run
python3 pipelines/backfill_glycoct.py --config pipelines/config.local.yaml
```

#### 多糖的完整结构从哪来

位移归属表只给出**每个残基的连接位点**（如 `→4)-β-D-Galp-(1→`），不含残基
之间的连接顺序。但多糖文献的作者已经用 2D NMR（HMBC/NOESY）把顺序推断好，
并把**结论**写在正文里，例如：

```
the connection of the main chain was
→[4)-β-D-Galp-(1]9→4,6)-β-D-Galp-(1→4)-α-D-GalpA-(1→...→4)-α/β-D-Glcp
```

因此 **ETL 不需要、也不应该重新解析 2D 谱** —— 只需抽出这句正文。解析流程：

1. `text_linear` 按栏线性化（双栏排版会把该句被另一栏切断）；
2. 定位结论句（"main chain/backbone ... was"），排除图注里的同名短语；
3. 展开重复块 `[X]n`，解析每个残基的位点、构型、环形式；
4. 同样提取支链定义（"branch chains were R1: … and R2: …"），按**糖基类型**
   匹配主链分支点后挂接 —— 匹配不唯一时退化为按顺序匹配，并在 `qc_notes`
   注明"待人工确认"，不静默臆造；
5. 交给 `glycoct.py` 生成完整序列编码 → `structure_level='complete'`。

实测效果（黄精多糖 SPR-1，J. Pharm. Anal. 2024）：

- 主链 15 残基（9×β-D-Galp + 4,6-β-D-Galp + 2×α-D-GalpA + α-D-Glcp +
  4,6-α-D-Glcp + 还原端 α-D-Glcp）
- 支链 R1 `β-D-Galp-(1→3)-β-D-Galp-(1→` → 挂到 4,6-β-D-Galp 的 O6；
  R2 `α-D-Glcp-(1→6)-α-D-Glcp-(1→` → 挂到 4,6-α-D-Glcp 的 O6
- 5 条链共 19 残基、18 条糖苷键，`composition = Gal12,Glc5,GalA2`，
  GlycoCT 通过 glypy 校验

表格归属的 12 个残基类型仍保留在 `residues` 表，用于与 152 个位移一一对应
（表格只列残基类型，正文连接式才给出聚合度）。若正文没有该结论句，
则如实退回 `composition_only`（如山楂多糖一文）。

---

## 5. 如何对接（下游平台 / 检索知识库 / AI 平台）

1. **直接查询**：任何 PostgreSQL 客户端按 §2.1 连接 `glycan_db`，从 8 张业务表读取结构化真值记录。
2. **向量检索（语义召回相似结构）**：库内含 embedding 列（pgvector），可用 SQL 近邻检索：
   ```sql
   SELECT s.iupac_short, 1 - (e.embedding <=> $1::vector) AS sim
   FROM sugars s JOIN nmr_experiments nx USING (sugar_id)
   JOIN nmr_shifts_1d e ON e.experiment_id = nx.experiment_id
   ORDER BY e.embedding <=> $1::vector
   LIMIT 10;
   ```
3. **结构比对（AI 推断 vs 真值）**：AI 平台给出候选结构的标准 GlycoCT 后可直接回查：
   ```sql
   -- 3.1 精确命中：候选结构是否已在真值库中（仅对 structure_level 为
   --     complete / repeat_unit 的记录有意义）
   SELECT s.sugar_id, s.iupac_short, s.structure_level, s.structure_confidence,
          s.first_seen_doi
   FROM sugars s
   WHERE s.glycoct = $1;

   -- 3.2 按组成筛选：连接顺序未知的杂多糖只能用组成检索
   SELECT s.iupac_short, s.composition, count(r.residue_id) AS n_residues
   FROM sugars s LEFT JOIN residues r ON r.sugar_id = s.sugar_id
   WHERE s.composition = $2 AND s.structure_level = 'composition_only'
   GROUP BY 1, 2;
   ```
   比对前先看 `structure_level`：`composition_only` 的记录没有结构编码，
   不能做结构级判定，只能提示"该文献报道过相同组成的多糖"。
4. **对接知识库**：将 `etl_reports/` 报告与结构化记录作为语料注入检索知识库，供上层 Agent/AI 平台引用。

---

## 6. 安全与规范

- 本仓库**禁止提交任何真实数据库口令 / API Key**：所有敏感值以占位符（`CHANGE_ME`、`sk-xxx`）或环境变量形式存在。
- 真实连接配置请复制为 `pipelines/config.local.yaml`（已被 `.gitignore` 排除）。
- Embedding 的 API Key 通过环境变量注入：`export SILICONFLOW_API_KEY=<你的key>`。
*（内容由AI生成，仅供参考）*
