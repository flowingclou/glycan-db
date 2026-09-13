# glycan-db 项目评估与优化建议

> 评估对象：`git@github.com:flowingclou/glycan-db.git`（本地 `main` = `origin/main` = `202bc3b`，8 次提交，无未提交改动）
> 评估方式：完整通读 `glycan_etl/*`、`pipelines/*`、`tests/*`、`db/*`；实跑自检与回归测试；在本地 PostgreSQL 上复现 `init.sql` / 迁移 / 写入幂等性；用 glypy 校验真实库存 GlycoCT。
> 结论：**项目流程清晰、设计取舍成熟、数据真实性纪律很好；当前主要风险不在"能不能解析"，而在"解析结果的身份（identity）不稳定"与"工程化护栏缺失"。**

---

## 一、项目流程（实际跑通的链路）

```
文献 SI PDF
   │
   ├─① 类型探测  batch_etl.peek_pdf_text → detect_pdf_type
   │     前 5 页打分（表格特征 vs 行内特征），无信号则深探到 12 页
   │     判为 table / inline / mixed；主路 0 条时自动补跑另一路（fallback）
   │
   ├─② 解析（两条独立路径）
   │     table_parser.py ── 词级坐标聚类重建行 → 分子量表(KV) + 位移归属表(残基×位置)
   │                        + 正文连接式（结论句）→ 支链提取与挂接 → 完整序列
   │     core.py ───────── 段落切分(按化合物标题) → 1D/2D 谱头与峰 → 物化/多糖性质
   │
   ├─③ 结构编码  glycoct.build_glycoct
   │     残基 + 连接 → 标准 GlycoCT（RES/LIN），并用 glypy 校验
   │     产出三级表达等级：complete / repeat_unit / composition_only
   │     **连不出顺序时拒绝伪造结构**，只给 composition
   │
   ├─④ 质控  Python 侧 R1/R2/R3 预检 + DB 侧触发器 R1-R10 + 全库校验函数
   │     flagged 默认仍入库并打标（可回溯），--skip-flagged 才丢弃
   │
   └─⑤ 入库  insert_record：literature UPSERT → sugars ON CONFLICT(glycoct)
         → 按 (结构, 文献) 先删后写派生行（experiments/shifts/2D/物化/多糖性质）
         
报告层：etl_reports/manifest_*.{csv,md}（批次汇总 + 0 条原因分类 + 结构等级）
向量层：embeddings.py 调 BGE-M3 回写 nmr_shifts_1d.embedding (vector 1024) + HNSW
```

**实测验证结果（不是读代码推测）**

| 验证项 | 结果 |
|---|---|
| `tests/run_tests.py` | 3/3 通过，回归测试 **31/31** 通过（README 写的是"11 项"，已过时） |
| `db/init.sql` 连续执行 2 次 | 均 `exit=0`，0 error，幂等性成立 |
| `init.sql` vs `v001..v005` 逐步执行 | 列集合、索引、触发器**完全一致**（无漂移，README 描述属实） |
| 真实库存 GlycoCT（SPR-1） | 19 残基 / 18 条键，glypy 校验 `True`；连接图是一棵无环连通树 ✓ |
| 真实库存 GlycoCT（Inulin） | glypy 校验 **False**（`KeyError: Could not translate fru`）—— 酮糖已知局限，仅记 warning |

---

## 二、值得肯定的设计（这些是有意为之的好取舍）

1. **不伪造结构**：`structure_level` 三级 + `composition_only` 时宁可留空 GlycoCT。这是真值库最要紧的纪律，很多同类项目在这里会直接编数据。
2. **占位编码带结构指纹 + DOI**：早期"空字符串撞 UNIQUE 导致结构互相覆盖"的坑已经踩过并修好，还配了防回退测试（`test_placeholder_unique_per_structure`）。README 甚至写明了"同一结构多文献会各占一条，待补全后归并"这个 trade-off。
3. **作者结论句优先于自行重解析 2D**：`→[4)-β-D-Galp-(1]9→…` 这类句子是作者用 HMBC/NOESY 推断好的结论，直接抽取比 ETL 自己重推 2D 谱可靠得多。这个判断非常专业。
4. **failed 数据保留而非丢弃**：QC flagged 默认入库打标，符合真值库"可回溯优先"的原则。
5. **词级坐标重建 + 双栏线性化 + x_tolerance=1.5**：解决了无框线对齐表格与"字间距紧密"排版两个真实痛点，且各自有注释解释为什么。
6. **回归测试针对"真实出现过的缺陷"**：测试文件头写明每个用例防的是哪个历史 bug，这是很成熟的工程习惯。

---

## 三、发现的问题（按优先级）

### P0-1 结构身份不稳定 → 解析器一升级就产生孤行（已实测复现）

`sugars.glycoct` 同时承担"结构表示"和"去重主键"两个角色，但它的取值会随解析路径变化：

```
同一份 PDF、同一份 config，跑两次 insert_record：
  第 1 次（解析器只到 composition_only）→ 写入 UNRESOLVED:poly:10.9999/audit:cd3341a03f0f529b
  第 2 次（解析器升级后给出真实序列）  → 写入 RES\n1b:a-dgal-HEX-1:5   ← 新行
  结果：sugars 7→8、nmr_experiments 8→9、residues 4→5（旧行永久留在库里）
```

`UNRESOLVED:...` 指纹由 `_structure_signature()` 拼接 `iupac_short / molecular_formula / 残基列表` 而成 —— 这三个字段**恰恰是最常随解析改进而变化的**。所以：

* README §2.4 承诺的"重复导入是安全的"只在解析器**逐字节稳定**时成立；
* 每次改进 parser 后重跑，都会给同一份文献留下新旧两条记录，谱图被劈成两半；
* 空的 `residues`（5 个 sugar 没有任何残基行）会让第二条记录的指纹漂移得更远。

**修复方向**：把"内容寻址"与"结构表示"拆开。

```sql
-- 1) glycoct 降级为普通列，唯一性改由内容键承担
ALTER TABLE sugars ADD COLUMN content_key TEXT;   -- 规范化的结构指纹（稳定哈希）
CREATE UNIQUE INDEX uq_sugars_content_key ON sugars(content_key);
ALTER TABLE sugars DROP CONSTRAINT sugars_glycoct_key;
-- 2) 解析器/版本变化留痕，而不是靠新行表达
ALTER TABLE sugars ADD COLUMN parser_version TEXT, ADD COLUMN last_parsed_at TIMESTAMPTZ;
-- 3) 每个结构被哪些文献观测到（多对多，替代"一结构一文献"的隐含假设）
CREATE TABLE sugar_sources (
  sugar_id BIGINT REFERENCES sugars ON DELETE CASCADE,
  source_id BIGINT REFERENCES literature,
  nmr_page INT, parser_version TEXT, structure_level TEXT, observed_at TIMESTAMPTZ DEFAULT now(),
  PRIMARY KEY (sugar_id, source_id)
);
-- 4) 同一 content_key 出现两种解析结果时不静默覆盖，落复核队列
ALTER TABLE sugars ADD COLUMN merge_status TEXT DEFAULT 'unique';  -- unique/conflict/merged
```

写入侧对应改为 `ON CONFLICT (content_key) DO UPDATE`，并在 `glycoct` 变化时写 `merge_status='conflict'` + 报告，而不是新增行。

### P0-2 `assignment_position` 存的是"残基码 + 位置"（已实测）

`table_parser.build_record()` 里 `Peak1D(assignment=f"{code} {posname}")` 同时填给了 `assignment` 和 `assignment_position`，线上真实数据：

```
GAt 2 | GE1,4 6a/6b | Rα 2 | GA1,2,4 5 | A1,5 4 | M 6a
```

后果三连：
1. 下游想按"H-1 / C-4 位点"筛选位移时要再做一次字符串解析；
2. `embeddings.py` 的 `build_text()` 会拼出 `"assigned to position GE1,4 6a/6b"` —— 语义文本质量下降，直接影响向量召回；
3. `qc_run_full_check()` 的 R4/R6 依赖 `assignment_position` 做位点语义判断，形同虚设。

**修复**：`assignment_position` 只存 `H-1` / `C-6a/6b` / `OCH3`，残基编码另开一列 `residue_code`（schema 加列 + 回填脚本），或至少在入库前 `split()` 一次。

### P1-1 解析粒度过粗导致的静默数据损失（已实测）

线上库：`α-D-Glcp` 只有 1 条 1H 位移，`α-D-Glcp-(1→4)-D-Glcp` 的 1H 有 1 条、2D 实验有 0 条相关峰。原因是 `core.parse_experiments()` 的段落上界：

```python
end = heads[idx + 1][0] if idx + 1 < len(heads) else pos + 1500   # ← 最后一个实验写死 1500 字符
```

最后一个 NMR 头之后没有下一个头时，只用 1500 字符；超出的峰被**静默丢弃**（无 warning、无计数）。同一段内如果下一个头出现得晚，13C 段落里的数字还会被 `RE_PEAK_1H` 按 ¹H 语义解析（数据串味）。

**修复**：`end = min(下一个头, 下一个标题/段落边界, 页面边界)`；并加"抽取峰数 vs 正则可匹配数值数"的对账，差集写入 `qc_notes` 计数而不是丢弃。

### P1-2 `pka` 字段永远为空（类型不匹配，静默）

`parse_physicochemical()` 对 `pKa` 做 `float()` 赋给 `pc.pka`，而 DDL 里 `physicochemical.pka` 是 `numeric(4,2)`；`pKa 3.5` 正好能过，但如 `pKa 12.3` 也在范围内 —— 真正的问题是 `pka` 从没被校验过是否真写了，而 QC 侧没有对应规则。建议加 DB 端 CHECK（0~14）与跨字段规则 R11（`pka` 非空则 `solubility` 应为酸性多糖描述等，或至少纳入 `qc_run_full_check`）。同类需核查：`melting_point_c` 区间、`optical_rotation` 符号与 D/L 的一致性。

### P1-3 `polysaccharide_props` 会累积重复行

`insert_record()` 删旧行的条件是 `WHERE sugar_id=%s AND source_id=%s`；同一结构被同一文献的两个段落命中时（`residues` 用 `ON CONFLICT (sugar_id, residue_seq)` 覆盖，但 `polysaccharide_props` 没有唯一约束、只会 INSERT），会产生多行。建议加 `UNIQUE (sugar_id, source_id, coalesce(repeat_unit_formula,''))` 类约束或按主键 UPSERT。

### P1-4 整个批次只在一个大事务里，重复导入会丢已入库的旧谱图

`insert_record()` 每个记录内部 `commit()`，但没有事务边界保护：任一条记录中途抛错（例如 `anomer` 触发 CHECK、`numeric` 溢出），该记录之前删除的派生行**不会恢复**，而 `batch_etl` 只 `log.error` 后继续下一份 PDF。表现是"重跑一次，某条记录的谱图反而变少了"。建议每条记录包一层 `SAVEPOINT`/独立事务，失败回滚该条并写入报告。

### P2 工程化与一致性

| 编号 | 问题 | 证据 | 建议 |
|---|---|---|---|
| P2-1 | 默认配置路径不可用，静默降级 dry-run | `--config` 默认 `etl_config.yaml`，仓库里只有 `pipelines/config.yaml`；在仓库根跑 `--dir X`（不带 `--config`）会打印"未找到配置…按 dry-run 执行"而不是报错 | 默认改为 `pipelines/config.yaml`，且路径可解析；找不到配置时 `exit 2` 而非静默降级 |
| P2-2 | 报告目录两套 | `REPORT_DIR="etl_reports"`（相对 CWD）→ 根目录 `etl_reports/` 有 54 个文件；`pipelines/etl_reports/` 又有一套 | 统一为 `--report-dir` 默认 `_REPO_ROOT/etl_reports`，跑完打印绝对路径 |
| P2-3 | `glycan_etl` 与 `pipelines/` 版本叙事不一致 | `__init__.py` = `3.2.0`；`docs/skill_guide.md` = "v1.0 基于 glycan_etl_v3.py v3.1"；`core.py` docstring = v3.1；`batch_etl.py` = v1.0 | 由 `__version__` 单一来源注入，文档引用变量 |
| P2-4 | 文档数字过时 | README 写"11 项回归测试"，实际 **31 项**；README 写"8 张业务表"，schema 实际 **9 张**（含 `spectrum_files`）；README §4.3 说 init.sql 是"四份迁移合并"，实际含 v005 | 把数字改为由测试输出/schema 生成，或至少在 CI 里校验 |
| P2-5 | 本地环境与文档不一致 | `pipelines/config.local.yaml` 指向 `glycan_final` / 用户 `hqx`；README 与 docker-compose 是 `glycan_db` / `glycan`。`embeddings.py` 更是**硬编码** `dbname="glycan_db"` 且不接受 `--config` | 环境差异写进 README（或直接用 docker-compose 作为唯一环境真相）；`embeddings.py` 改为读同一 config |
| P2-6 | 向量能力实际为 0 | 实测 `nmr_shifts_1d`：`count(*)=329`，`count(embedding)=0` | 跑一次 `embeddings.py` 并加"embedding 覆盖率"到报告；README §5 的 SQL 目前必然返回空集 |
| P2-7 | 无 CI | 无 `.github/workflows`；回归测试全靠人工 `python3 tests/run_tests.py` | 加 GitHub Actions：`ruff` + `run_tests.py` + 用 service container 跑 `init.sql` 两遍 + 约束测试 |
| P2-8 | 无 `pyproject.toml` / 无打包元数据 | 只有 `requirements.txt`（无锁定版本、无 hash） | 加 `pyproject.toml`（`[project]` + 可选 `[project.optional-dependencies] dev`），用 `uv`/`pip-tools` 生成锁文件 |
| P2-9 | 硬编码坐标常量（移植性主要瓶颈） | `table_parser.py`：残基列 `< 90`、IUPAC 列 `88..195`、数值列 `>= 220`、KV 标签列 `< 160` / 值列 `160..330`、`hc_anchor` 兜底 `205.0` | 改为按表头词坐标推导列锚点（已有 `col_anchor` 机制，扩展成通用版）+ 按"相对列宽比例"；每份新 PDF 先跑一次"列检测自检"，失败的版面进复核队列 |
| P2-10 | 一处 `and`/`or` 优先级 bug | `table_parser.py:414` `re.search(r"6a/6b", joined) or re.search(r"-OMe", joined) and re.search(r"\b5\b", joined)` → 只有"6a/6b"或"带 5 的 -OMe"才被认作位置表头；实测 `1 2 3 4 6 -OMe` 返回 False（本意应为 True）。第 428-430 行的"仅数字标题也可"兜底分支因此**永不可达**（425 行已提前 return） | 补括号；补一个"无 6a/6b、有 -OMe"的版式测试 |
| P2-11 | 无交叉一致性校验（静默错误的主要来源） | 单糖 `β-D-GlcpNAc` 的 `molecular_formula` 来自 `MONO_FORMULA` 的 **Glc 值** `C6H12O6`（N-乙酰氨基糖应有 N） | `MONO_FORMULA` 补全 `GalA/Rha/Ara/Xyl/GlcN/GalN/MurNAc` 等；加"由 GlycoCT 推导分子式 ↔ 文献分子式"对账，超差即 `flagged` |
| P2-12 | 未使用的导入 / 死代码 | `core.py`：`sys`；`batch_etl.py`：`re, sys, csv` 中 `re/sys` 未用；`table_parser.py`：`Counter, Tuple, defaultdict, Any`；`table_parser.cluster_rows` 里 `if best is None or best["top"] == 0 and rows: pass` 是空语句 | 引入 `ruff` 并开启 `F401/F841` |

---

## 四、优化建议（按顺序执行）

### 第 1 步（本周，直接决定数据可信度）
1. **拆分"结构表示"与"去重主键"**（P0-1）：加 `content_key` + `sugar_sources` + `merge_status`，`ON CONFLICT (content_key)`；写一份 `v006_identity.sql` 迁移 + 回填脚本（现有 8 条记录可人工确认后归并）。
2. **修 `assignment_position` 语义**（P0-2）：加 `residue_code` 列，回填 `UPDATE ... SET assignment_position = split_part(...)`，重跑 `embeddings.py`。
3. **`--config` 默认值 + 缺配置即报错**（P2-1），并让报告目录固定（P2-2）。这一条能直接消除"以为入库了其实只是 dry-run"的隐性事故。
4. **跑一次向量化并核对维度**（P2-6），把 embedding 覆盖率写进 manifest。

### 第 2 步（两周，消除静默错误）
5. `parse_experiments` 段落边界改为真实边界 + 抽取对账（P1-1）。
6. 交叉校验层：GlycoCT 推导分子式 / 残基组成 vs 文献值；`pka`、熔点、旋光度区间规则（P1-2、P2-11）。
7. 事务与幂等：按记录 `SAVEPOINT`，`polysaccharide_props` 加唯一约束（P1-3、P1-4）。
8. 补 `tests/` ：`pytest` 化 + SQL 约束/迁移测试（`init.sql` 跑两遍 + 触发器断言，我已验证这套做法可行）+ 真实 PDF 的 golden JSON 快照。

### 第 3 步（一个月，扩展性与可维护性）
9. **置信度分级**：`structure_confidence` 目前只有 `confirmed_2d/confirmed_1d/reported`，建议扩展为 `curated > parsed_sequence > table_assigned > composition_only`，并写清"下游哪一级可以直接按结构比对"。
10. **列检测通用化**（P2-9 + P2-10）：让 table_parser 先输出"本次识别到的列锚点与置信度"，人工/程序确认后再解析；这是把"8 篇论文可用"变成"可持续加文献"的关键一步。
11. **CI + 打包 + 版本单一来源**（P2-3、P2-7、P2-8）。
12. **文档与代码对齐**（P2-4、P2-5），并把本文档中的实测数字写进 README 的"当前状态"小节。

---

## 五、一句话总结

引擎的**判断力**（不伪造结构、优先采信作者结论句、失败数据保留）已经超过多数同类项目；当前短板是**工程护栏**——结构身份靠可变字符串寻址、派生数据靠"先删后写"却没有事务保护、配置默认值会让入库静默降级为 dry-run、向量能力是空的。把 P0-1 / P0-2 / P1-1 / P1-4 四件事做掉，这套 ETL 就从"能跑通 8 篇文献"升级为"可以持续加文献、且每次都敢信结果"。
