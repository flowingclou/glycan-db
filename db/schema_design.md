---
AIGC:
    Label: "1"
    ContentProducer: 001191440300708461136T1XGW3
    ProduceID: 463222f6ea5b7ebcc547bc2191331f46_02ced6baae7511f1b128525400f8a581
    ReservedCode1: sDuQiKgReROiLnxlY1riUjF4JGnpWGb8jBqwlLjPqB9L52umywhdaFg2DbkG711n6pHQWFZx9MfVjZDz/Svh1aJTPXknUmwZ0Vl//3FtZK/QtODKatkesM6APELtKlz26VBbEowA2NT7Cb55a5togEbz+q+0k6gWKGbewKfL3/Lcqd6T6LKWSKJg+vk=
    ContentPropagator: 001191440300708461136T1XGW3
    PropagateID: 463222f6ea5b7ebcc547bc2191331f46_02ced6baae7511f1b128525400f8a581
    ReservedCode2: sDuQiKgReROiLnxlY1riUjF4JGnpWGb8jBqwlLjPqB9L52umywhdaFg2DbkG711n6pHQWFZx9MfVjZDz/Svh1aJTPXknUmwZ0Vl//3FtZK/QtODKatkesM6APELtKlz26VBbEowA2NT7Cb55a5togEbz+q+0k6gWKGbewKfL3/Lcqd6T6LKWSKJg+vk=
---

# 糖类数据库（单糖 / 寡糖 / 多糖）Schema 设计方案

> 参考 NMRexp（Scientific Data, 2025）小分子 NMR 数据库构建范式，针对糖类结构表征特点（糖苷键连接、异头构型、2D NMR 归属、多糖宽峰）做的专项改造设计。
> 版本：v1.0 · 设计日期：2026-09-12

---

## 1. 总体架构

### 1.1 设计原则

1. **分层建模**：单糖 / 寡糖 / 多糖结构粒度和谱图复杂度差异大，拆分为独立主表，不硬塞一个 schema。
2. **结构主键用糖链编码**（GlycoCT / WURCS），SMILES 仅作可选交叉字段。
3. **"结构-谱图对" + 完整溯源**：每条谱图记录可回溯到源文献 DOI 与页码。
4. **宁缺毋滥的质控**：全部清洗规则代码化，每条记录带确认等级标记。
5. **原始谱图文件与解析结果分离存储**：既存结构化峰列表，也保留原始 1D/2D 谱文件。

### 1.2 架构总览（数据流）

```
文献 PDF/SI
   │  ① 版面检测 + 文本抽取
   ▼
结构抽取 ──→ 糖链编码(GlycoCT/WURCS/IUPAC) + 残基表
   │  ② 规则 + LLM 混合解析
   ▼
NMR 文本解析 ──→ 1D峰表 + 2D相关峰表
   │  ③ 糖专属质控规则校验
   ▼
标准化入库 ──→ 糖主表 / 残基表 / 谱图表 / 物化性质表 / 多糖专属表
   │
   ▼
Parquet 导出（供 AI 训练）+ PostgreSQL 在线查询
```

---

## 2. 糖链编码方案（核心）

### 2.1 编码体系选择

| 编码 | 用途 | 定位 |
|---|---|---|
| **GlycoCT**（XML） | 结构主键 / 规范化存储 | 能表达分支、连接位置、异头构型、环形开链，是糖类事实标准，推荐为主键 |
| **WURCS** | 主键备选 / 与 GlyTouCan 互操作 | GlyTouCan（糖类结构登记库）官方编码，可做跨库关联 |
| **IUPAC 缩写字序列** | 人类可读字段 | 如 `Glcα1-4Glc`、`GalNAcβ1-4(Neu5Acα2-3)Galβ1-4Glc`，用于展示与检索 |
| **SMILES / InChIKey** | 可选交叉字段 | 糖用 SMILES 存在开环/构型歧义，仅作交叉参考，不参与主键 |

### 2.2 关键约束

- **主键生成规则**：`GlycoCT 规范化字符串`（去空白、统一大小写）的 SHA-256 哈希 或直接存原文，建议直接存规范原文 + 哈希双字段，哈希用于去重与索引。
- **残基粒度**：每条糖结构拆成 `residues` 子表，记录每个残基的连接位置（parent_carbon → child_anomeric_carbon）、异头构型（α/β）、环形式（pyranose/furanose）。
- **结构唯一性**：同一 GlycoCT 唯一对应一条 `sugars` 主记录，不同文献报到同一结构时共用主键、挂多条谱图记录。
- **多糖特例**：多糖存"重复单元"的结构编码（GlycoCT 支持重复单元语法 `B|b:...` 标记），聚合度/分子量放 `polysaccharide_props` 表。

---

## 3. 表结构设计（9 张表）

### 3.1 `sugars` —— 糖分子主表（核心）

| 字段 | 类型 | 约束 | 说明 |
|---|---|---|---|
| sugar_id | BIGSERIAL | PK | 内部自增 ID |
| sugar_type | TEXT | NOT NULL, CHECK IN ('mono','oligo','poly') | 单糖/寡糖/多糖 |
| glycoct | TEXT | NOT NULL, UNIQUE | 规范化 GlycoCT 编码（主键依据） |
| glycoct_hash | CHAR(64) | NOT NULL, UNIQUE | SHA-256，用于去重索引 |
| wurcs | TEXT | | WURCS 编码（GlyTouCan 互操作） |
| iupac_short | TEXT | NOT NULL | IUPAC 缩写字序列，如 Glcα1-4Glc |
| smiles | TEXT | | 可选交叉编码 |
| inchikey | TEXT | | 可选交叉编码 |
| molecular_formula | TEXT | | 分子式 |
| molecular_weight | NUMERIC(10,4) | | 分子量（Da） |
| glytoucan_id | TEXT | | GlyTouCan 注册 ID |
| monosaccharidedb_id | TEXT | | MonosaccharideDB ID |
| structure_confidence | TEXT | NOT NULL DEFAULT 'reported', CHECK IN ('confirmed_2d','confirmed_1d','inferred','reported') | 结构确认等级 |
| stereochemistry_defined | BOOLEAN | NOT NULL DEFAULT FALSE | 立体化学是否完全明确（α/β、D/L） |
| first_seen_doi | TEXT | | 首次出现的源文献 DOI |
| created_at / updated_at | TIMESTAMPTZ | NOT NULL | 时间戳 |

### 3.2 `residues` —— 糖残基表（结构组成）

| 字段 | 类型 | 约束 | 说明 |
|---|---|---|---|
| residue_id | BIGSERIAL | PK | |
| sugar_id | BIGINT | FK → sugars.sugar_id, NOT NULL | 所属糖 |
| residue_seq | INT | NOT NULL | 残基序号（按 IUPAC 从还原端编号） |
| monosaccharide_name | TEXT | NOT NULL | 单糖名，如 Glc、GalNAc、Neu5Ac |
| ring_form | TEXT | CHECK IN ('p','f','open','unknown') | 吡喃/呋喃/开链 |
| anomer | TEXT | CHECK IN ('a','b','unknown') | 异头构型 |
| is_reducing_end | BOOLEAN | DEFAULT FALSE | 是否还原端 |
| parent_carbon | INT | | 糖苷键母体碳位（如 4） |
| linkage_branch | INT | | 分支点（0=线性） |

### 3.3 `nmr_experiments` —— 谱图实验记录（谱图头部）

| 字段 | 类型 | 约束 | 说明 |
|---|---|---|---|
| experiment_id | BIGSERIAL | PK | |
| sugar_id | BIGINT | FK → sugars, NOT NULL | 关联结构 |
| source_id | BIGINT | FK → literature, NOT NULL | 溯源 |
| nmr_type | TEXT | NOT NULL, CHECK IN ('1H','13C','2D','19F','31P','other') | 谱类型 |
| nucleus_list | TEXT[] | | 所含核素 |
| solvent | TEXT | NOT NULL | 溶剂，如 CDCl3、D2O |
| frequency_mhz | NUMERIC(7,2) | | 仪器 ¹H 频率（MHz） |
| temperature_c | NUMERIC(6,2) | | 温度 ℃ |
| ph | NUMERIC(4,2) | | pH（水相谱图必填） |
| assignment_level | TEXT | CHECK IN ('full','partial','residue_level','not_assigned') | 归属粒度（多糖常为 residue_level） |
| spectrum_file_id | BIGINT | FK → spectrum_files | 原始谱文件 |
| qc_status | TEXT | NOT NULL DEFAULT 'pending', CHECK IN ('pending','passed','flagged','rejected') | 质控状态 |

### 3.4 `nmr_shifts_1d` —— 一维化学位移峰列表

| 字段 | 类型 | 约束 | 说明 |
|---|---|---|---|
| shift_id | BIGSERIAL | PK | |
| experiment_id | BIGINT | FK → nmr_experiments, NOT NULL | 所属实验 |
| nucleus | TEXT | NOT NULL, CHECK IN ('1H','13C') | 核素 |
| shift_ppm | NUMERIC(7,3) | NOT NULL | 化学位移 |
| multiplicity | TEXT | | 多重峰类型（s/d/t/q/m/br 等） |
| j_coupling_hz | NUMERIC(7,2) | | J 耦合常数 |
| integration | NUMERIC(6,2) | | 积分/氢计数（¹H） |
| assignment_residue | TEXT | | 归属残基（如 Glc-A） |
| assignment_position | TEXT | | 归属位置（如 H1、C1、C4 等） |
| assignment_confidence | TEXT | CHECK IN ('confirmed_2d','predicted','unassigned') | 归属置信度 |
| is_anomeric | BOOLEAN | DEFAULT FALSE | 是否异头碳/氢（糖专属关键标记） |

### 3.5 `nmr_correlations_2d` —— 二维谱相关峰表（寡/多糖核心扩展）

| 字段 | 类型 | 约束 | 说明 |
|---|---|---|---|
| corr_id | BIGSERIAL | PK | |
| experiment_id | BIGINT | FK → nmr_experiments, NOT NULL | |
| experiment_2d | TEXT | NOT NULL, CHECK IN ('HSQC','HMBC','COSY','TOCSY','NOESY','ROESY') | 二维谱类型 |
| proton_shift_ppm | NUMERIC(7,3) | NOT NULL | ¹H 位移 |
| hetero_shift_ppm | NUMERIC(7,3) | | ¹³C 位移（HSQC/HMBC） |
| atom_pair | TEXT | NOT NULL | 相关原子对，如 H1-C1、H1-C4' |
| residue_from | TEXT | | 起始残基 |
| residue_to | TEXT | | 目标残基（用于 HMBC/NOESY 糖苷键判定） |
| linkage_evidence | TEXT | | 该相关峰是否为糖苷键连接证据（如 HMBC H1→C4'） |
| intensity | NUMERIC(6,2) | | 相对强度（NOESY/ROESY） |

### 3.6 `physicochemical` —— 物化性质表

| 字段 | 类型 | 约束 | 说明 |
|---|---|---|---|
| prop_id | BIGSERIAL | PK | |
| sugar_id | BIGINT | FK → sugars, NOT NULL | |
| optical_rotation | NUMERIC(8,3) | | 旋光度 [α]D |
| optical_rotation_condition | TEXT | | 旋光度测试条件（溶剂/浓度/温度） |
| melting_point_c | NUMERIC(7,2) | | 熔点 |
| solubility | TEXT | | 溶解性 |
| pka | NUMERIC(5,2) | | 酸性糖 pKa（如糖醛酸/唾液酸） |
| source_id | BIGINT | FK → literature | 溯源 |

### 3.7 `polysaccharide_props` —— 多糖专属性质表

| 字段 | 类型 | 约束 | 说明 |
|---|---|---|---|
| poly_id | BIGSERIAL | PK | |
| sugar_id | BIGINT | FK → sugars, NOT NULL | |
| repeat_unit_formula | TEXT | | 重复单元分子式 |
| degree_of_polymerization | NUMERIC(10,2) | | 平均聚合度 |
| molecular_weight_mn | NUMERIC(12,2) | | 数均分子量 |
| molecular_weight_mw | NUMERIC(12,2) | | 重均分子量 |
| polydispersity | NUMERIC(6,3) | | 分散度 PDI |
| monosaccharide_ratio | JSONB | | 残基摩尔比，如 {"Glc":0.7,"Gal":0.3}（来自甲基化/GC-MS） |
| backbone | TEXT | | 主链描述 |
| branching | TEXT | | 分支结构描述 |

### 3.8 `literature` —— 文献溯源表

| 字段 | 类型 | 约束 | 说明 |
|---|---|---|---|
| source_id | BIGSERIAL | PK | |
| doi | TEXT | NOT NULL, UNIQUE | 源文献 DOI |
| journal | TEXT | | 期刊 |
| year | INT | | 年份 |
| title | TEXT | | 标题 |
| authors | TEXT | | 作者 |
| supporting_info | BOOLEAN | DEFAULT FALSE | 是否来自 SI 补充材料 |
| nmr_page | INT | | NMR 数据所在页 |
| structure_page | INT | | 结构所在页 |

### 3.9 `spectrum_files` —— 原始谱图文件表

| 字段 | 类型 | 约束 | 说明 |
|---|---|---|---|
| file_id | BIGSERIAL | PK | |
| experiment_id | BIGINT | FK → nmr_experiments | 关联实验（可空，允许多谱一文件） |
| file_format | TEXT | CHECK IN ('jcamp','nmredata','mnova','fid','pdf','png','other') | 文件格式 |
| file_path | TEXT | NOT NULL | 存储路径 |
| file_hash | CHAR(64) | NOT NULL | SHA-256 去重 |
| nmredata_valid | BOOLEAN | DEFAULT FALSE | 是否符合 NMReDATA 标准 |
| license | TEXT | | 版权/再分发许可 |

---

## 4. PostgreSQL DDL

```sql
-- ============ 糖类数据库 DDL（PostgreSQL 16+） ============
-- 扩展：JSONB 已内置；如需化学结构相似度可后续加 RDKit 扩展

-- 4.1 糖分子主表
CREATE TABLE sugars (
    sugar_id            BIGSERIAL PRIMARY KEY,
    sugar_type          TEXT NOT NULL CHECK (sugar_type IN ('mono','oligo','poly')),
    glycoct             TEXT NOT NULL UNIQUE,
    glycoct_hash        CHAR(64) NOT NULL UNIQUE,
    wurcs               TEXT,
    iupac_short         TEXT NOT NULL,
    smiles              TEXT,
    inchikey            TEXT,
    molecular_formula   TEXT,
    molecular_weight    NUMERIC(10,4),
    glytoucan_id        TEXT,
    monosaccharidedb_id TEXT,
    structure_confidence TEXT NOT NULL DEFAULT 'reported'
                        CHECK (structure_confidence IN ('confirmed_2d','confirmed_1d','inferred','reported')),
    stereochemistry_defined BOOLEAN NOT NULL DEFAULT FALSE,
    first_seen_doi      TEXT,
    created_at          TIMESTAMPTZ NOT NULL DEFAULT now(),
    updated_at          TIMESTAMPTZ NOT NULL DEFAULT now()
);
CREATE INDEX idx_sugars_type ON sugars(sugar_type);
CREATE INDEX idx_sugars_iupac ON sugars(iupac_short);

-- 4.2 糖残基表
CREATE TABLE residues (
    residue_id     BIGSERIAL PRIMARY KEY,
    sugar_id       BIGINT NOT NULL REFERENCES sugars(sugar_id) ON DELETE CASCADE,
    residue_seq    INT NOT NULL,
    monosaccharide_name TEXT NOT NULL,
    ring_form      TEXT CHECK (ring_form IN ('p','f','open','unknown')),
    anomer         TEXT CHECK (anomer IN ('a','b','unknown')),
    is_reducing_end BOOLEAN DEFAULT FALSE,
    parent_carbon  INT,
    linkage_branch INT DEFAULT 0,
    UNIQUE (sugar_id, residue_seq)
);

-- 4.3 谱图实验记录
CREATE TABLE nmr_experiments (
    experiment_id    BIGSERIAL PRIMARY KEY,
    sugar_id         BIGINT NOT NULL REFERENCES sugars(sugar_id) ON DELETE CASCADE,
    source_id        BIGINT NOT NULL REFERENCES literature(source_id),
    nmr_type         TEXT NOT NULL CHECK (nmr_type IN ('1H','13C','2D','19F','31P','other')),
    nucleus_list     TEXT[],
    solvent          TEXT NOT NULL,
    frequency_mhz    NUMERIC(7,2),
    temperature_c    NUMERIC(6,2),
    ph               NUMERIC(4,2),
    assignment_level TEXT NOT NULL DEFAULT 'partial'
                     CHECK (assignment_level IN ('full','partial','residue_level','not_assigned')),
    spectrum_file_id BIGINT REFERENCES spectrum_files(file_id),
    qc_status        TEXT NOT NULL DEFAULT 'pending'
                     CHECK (qc_status IN ('pending','passed','flagged','rejected'))
);
CREATE INDEX idx_exp_sugar ON nmr_experiments(sugar_id);
CREATE INDEX idx_exp_solvent ON nmr_experiments(solvent);

-- 4.4 一维化学位移峰
CREATE TABLE nmr_shifts_1d (
    shift_id            BIGSERIAL PRIMARY KEY,
    experiment_id       BIGINT NOT NULL REFERENCES nmr_experiments(experiment_id) ON DELETE CASCADE,
    nucleus             TEXT NOT NULL CHECK (nucleus IN ('1H','13C')),
    shift_ppm           NUMERIC(7,3) NOT NULL,
    multiplicity        TEXT,
    j_coupling_hz       NUMERIC(7,2),
    integration         NUMERIC(6,2),
    assignment_residue  TEXT,
    assignment_position TEXT,
    assignment_confidence TEXT CHECK (assignment_confidence IN ('confirmed_2d','predicted','unassigned')),
    is_anomeric         BOOLEAN DEFAULT FALSE
);
CREATE INDEX idx_shift_exp ON nmr_shifts_1d(experiment_id);

-- 4.5 二维谱相关峰
CREATE TABLE nmr_correlations_2d (
    corr_id          BIGSERIAL PRIMARY KEY,
    experiment_id    BIGINT NOT NULL REFERENCES nmr_experiments(experiment_id) ON DELETE CASCADE,
    experiment_2d    TEXT NOT NULL CHECK (experiment_2d IN ('HSQC','HMBC','COSY','TOCSY','NOESY','ROESY')),
    proton_shift_ppm NUMERIC(7,3) NOT NULL,
    hetero_shift_ppm NUMERIC(7,3),
    atom_pair        TEXT NOT NULL,
    residue_from     TEXT,
    residue_to       TEXT,
    linkage_evidence BOOLEAN DEFAULT FALSE,
    intensity        NUMERIC(6,2)
);
CREATE INDEX idx_2d_exp ON nmr_correlations_2d(experiment_id);

-- 4.6 物化性质
CREATE TABLE physicochemical (
    prop_id          BIGSERIAL PRIMARY KEY,
    sugar_id         BIGINT NOT NULL REFERENCES sugars(sugar_id) ON DELETE CASCADE,
    optical_rotation NUMERIC(8,3),
    optical_rotation_condition TEXT,
    melting_point_c  NUMERIC(7,2),
    solubility       TEXT,
    pka              NUMERIC(5,2),
    source_id        BIGINT REFERENCES literature(source_id)
);

-- 4.7 多糖专属性质
CREATE TABLE polysaccharide_props (
    poly_id               BIGSERIAL PRIMARY KEY,
    sugar_id              BIGINT NOT NULL REFERENCES sugars(sugar_id) ON DELETE CASCADE,
    repeat_unit_formula   TEXT,
    degree_of_polymerization NUMERIC(10,2),
    molecular_weight_mn   NUMERIC(12,2),
    molecular_weight_mw   NUMERIC(12,2),
    polydispersity        NUMERIC(6,3),
    monosaccharide_ratio  JSONB,
    backbone              TEXT,
    branching             TEXT
);

-- 4.8 文献溯源
CREATE TABLE literature (
    source_id       BIGSERIAL PRIMARY KEY,
    doi             TEXT NOT NULL UNIQUE,
    journal         TEXT,
    year            INT,
    title           TEXT,
    authors         TEXT,
    supporting_info BOOLEAN DEFAULT FALSE,
    nmr_page        INT,
    structure_page  INT
);

-- 4.9 原始谱图文件
CREATE TABLE spectrum_files (
    file_id         BIGSERIAL PRIMARY KEY,
    experiment_id   BIGINT REFERENCES nmr_experiments(experiment_id),
    file_format     TEXT CHECK (file_format IN ('jcamp','nmredata','mnova','fid','pdf','png','other')),
    file_path       TEXT NOT NULL,
    file_hash       CHAR(64) NOT NULL UNIQUE,
    nmredata_valid  BOOLEAN DEFAULT FALSE,
    license         TEXT
);
```

---

## 5. 糖专属质控规则（清洗逻辑，代码化）

| # | 规则 | 说明 |
|---|---|---|
| R1 | **异头碳/氢位移区间校验** | ¹³C 异头碳 90–110 ppm（醛糖）、¹H 异头氢 4.2–5.8 ppm，越界标 `flagged` |
| R2 | **J 耦合常数合理性** | ¹J(C-H) 150–180 Hz；³J(H1-H2) α-葡萄糖 ~3.5 Hz、β-葡萄糖 ~7.5 Hz，用于交叉验证构型 |
| R3 | **构型-耦合一致性** | 若 anomer=β 但 J(H1,H2)>6.5 Hz 区间不符，标记冲突 |
| R4 | **OH 信号处理** | 可交换质子（OH）在 D2O 中应缺失，在 DMSO-d6 中出现；按溶剂差异化校验 |
| R5 | **糖苷键一致性** | HMBC H1→C4' 相关峰应与 residues 表中连接位置一致，不一致标 `flagged` |
| R6 | **位移单调性** | 同残基内归属位移序列非单调（排除正常重叠）则标记 |
| R7 | **氢计数守恒** | ¹H 积分总数 vs 分子式氢数（扣除交换氢）一致性校验 |
| R8 | **多糖残基级归属** | 多糖仅允许 `residue_level` 归属，禁止伪造原子级归属；缺失即 `not_assigned` |
| R9 | **溯源强制** | 任何谱图记录缺少 DOI 直接 `rejected`（宁缺毋滥） |
| R10 | **立体化学子集** | `stereochemistry_defined=true` 提供高可信过滤子集（复刻 NMRexp 子集设计） |

---

## 6. 文件格式与存储方案

| 用途 | 格式 | 说明 |
|---|---|---|
| 在线查询/关系管理 | PostgreSQL（如上 DDL） | 主存储 |
| AI 训练数据集 | Parquet（列式） | 从 PostgreSQL 导出宽表，每行一个"结构-谱图对" |
| 原始谱图 | JCAMP-DX（优先）/ NMReDATA | 开放标准，建议推进 NMReDATA 标准化 |
| 结构交换 | GlycoCT XML / WURCS | 跨库互操作 |
| 发布 | Zenodo（数据集 + 代码，含清洗规则） | 复刻 NMRexp 的可复现发布模式 |

---

## 7. 实施路线建议

1. **MVP 阶段**：先做 `sugars` + `nmr_experiments` + `nmr_shifts_1d` + `literature` 四表，以单糖为主，锁定 1–2 个期刊时间窗验证流水线。
2. **扩展阶段**：加入 2D 相关峰表与多糖专属表，引入 GlycoCT 规范化管线。
3. **质控阶段**：实现 R1–R10 规则集，做 300 条人工抽样验证（复刻 NMRexp 的精度评估方法）。
4. **发布阶段**：导出 Parquet + Zenodo 发布 + 清洗代码开源。
*（内容由AI生成，仅供参考）*
