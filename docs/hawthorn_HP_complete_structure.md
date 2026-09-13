# 山楂多糖 HP 的「完整结构」转录（Fig. 3E，图片格式）

> 来源：Zhang Y. et al., *Structural characterisation of hawthorn polysaccharide and its mechanisms of action against non-alcoholic fatty liver disease*, **Int. J. Biol. Macromol. 319 (2025) 145713**, DOI `10.1016/j.ijbiomac.2025.145713`
> 转录依据：Fig. 3 面板 **(E) Predicted structure of HP**（该面板整体为**位图**，PDF 文本层中不存在）
> 证据图：`docs/evidence/HP_fig3E_HG.png`、`HP_fig3E_RGI_arabinan.png`、`HP_fig3_full.png`
> 转录状态：**已人工核对**（600/900/1000/1600/1900/2400 dpi 多次渲染交叉确认）

---

## 1. 转录结果（原文写法）

**① HG 域（同型半乳糖醛酸聚糖主链）**

```
HG:  ···→4)-α-D-GalpA-(1→4)-α-D-GalpA-(1-[→4)-α-D-GalpAR-(1-]80→4)-α-D-GalpA-(1→···
     R = H or methyl ester
     Galacturonic acid has a small amount of branched chain O-2 and O-3
```

**② RG-I 域（鼠李糖半乳糖醛酸聚糖 I 主链 + 半乳聚糖侧链 R₁）**

```
RG-I: ···-[-→4)-α-D-GalpA-(1→2)-α-L-Rhap-(1-]2→4)-α-D-GalpA-(1→2)-α-L-Rhap-(1→···
                                              │4
                                              ↑
                                             R₁

R₁ = ···→3)-β-D-Galp-(1→3)-β-D-Galp-(1→···
              │6
              ↑1
           β-D-Galp
              │6
              ↑
              ⋮
```

**③ 阿拉伯聚糖链（arabinan）**

```
···→5)-α-L-Araf-(1→5)-α-L-Araf-(1→···
                    │3
                    ↑
                    ⋮
```

---

## 2. 转写为数据库可用的结构表示

### 2.1 域级组成

| 域 | 主链单元 | 连接方式 | 侧链/分支 | 原文给出的聚合度 |
|---|---|---|---|---|
| HG | α-D-GalpA | (1→4) | R = H 或甲酯；少量 O-2 / O-3 分支 | 括号下标 **80** |
| RG-I | [-→4)-α-D-GalpA-(1→2)-α-L-Rhap-(1-] | 二糖重复单元 | Rha 的 **O4** 接 R₁ | 下标 **2**（示重复） |
| R₁（半乳聚糖） | β-D-Galp | (1→3) | 每 6 位再接 β-D-Galp，可继续延伸 | ⋮ 省略 |
| 阿拉伯聚糖 | α-L-Araf | (1→5) | 3 位分支，可继续延伸 | ⋮ 省略 |

### 2.2 与 GlycoCT 的对应关系（供 `glycan_etl/glycoct.py` 使用）

```text
RES
1b:a-dgal-HEX-1:5|6:a          # α-D-GalpA（主链重复单元）
2b:a-lrha-HEX-1:5              # α-L-Rhap
3b:b-dgal-HEX-1:5              # β-D-Galp（R₁ 链）
4b:a-lara-PEN-1:4              # α-L-Araf
LIN
1:1o(4+1)1d                    # HG: GalA-(1→4)-GalA 重复
2:1o(2+1)2d                    # RG-I: Rha-(1→2)-GalA
3:2o(4+1)1d                    # RG-I: GalA-(1→4)-Rha
4:2o(4+1)3d                    # R₁ 经 Rha O4 接入
5:3o(3+1)3d                    # R₁ 内部 (1→3)
6:4o(5+1)4d                    # 阿拉伯聚糖 (1→5)
```

> 注意：以上仅为**重复单元/域级**编码，**不是**逐残基的完整序列——
> 原图本身就是"可能结构要素"（possible structural elements）的示意，
> 作者亦未给出残基到残基的唯一连接顺序。落库时必须标注为域级/示意级，
> 不可标记为 `complete`。

---

## 3. 为什么现有 ETL 管线拿不到这段结构（三个独立原因）

### 原因 1：Fig. 3E 整体是**位图**，文本层里不存在

实测：对第 6 页做 `get_text("dict")`，**只有图注一个文本块**，
面板 (E) 内的全部文字（`HG:`、`RG-I:`、`R₁=`、`R=H or methyl ester` …）
均**无文本行**，扫描 y>470 区域得到 0 条文本行。
→ `page.extract_text()` 无论用什么 `x_tolerance` 都取不到，纯文本管线必然漏掉。

### 原因 2：正文其实**写了**结论，但措辞不匹配触发词

论文 3.2 节原文（**文本层里存在**）：

> "These results suggest that HP is mainly composed of **a large number of HG domains
> and a small number of RG-I domains with side chains**. The possible structural
> elements are shown in **Fig. 3E**."

而 `table_parser.extract_linkage_sequence()` 的触发词只有：

```python
r"(?:main\s+chain|backbone)[^.]{0,40}?\b(?:was|were|is)\b\s*([→\[][^.]{10,800})"
```

—— 要求 `main chain/backbone` + `→` 或 `[` 开头。这句是 `composed of ... HG domains`，
**不含箭头**，因此不触发。这是**纯文本就能拿到**的结论句，漏掉纯属性问题，不是信息缺失。

### 原因 3：Table 2 的 ¹³C 行虽被解析到，但结构等级仍无法判定

实测 `parse_shift_table` 正确抓到 **16 个残基、82 个 ¹H 峰、88 个 ¹³C 峰**
（含 170.66/170.67/170.68/175.01/175.02/175.4 羧基碳，98–110 ppm 异头碳 12 个），
数据本身没丢。但 Table 2 **只有"每个残基被取代的位点"**，
没有残基 A 连到残基 B 的顺序信息，所以 `build_glycoct` 正确地退回 `composition_only`。

---

## 4. 论文内部数据不一致（人工复核时需注意）

| 数据源 | 内容 | 冲突点 |
|---|---|---|
| 单糖组成（Table 1 / 正文 §3.6） | GalA 85.101 %、**Man 4.506 %**、Gal 3.983 % | 合计仅 **93.6 %**；**Man 未出现在任何残基表中**；**无 Ara / Rha / GlcA 的百分比** |
| 甲基化分析（Fig. 1E / §3.1） | 4-GalA(p) 为主，另有 t-GalA(p)、4-Glc(p)、2,4-GalA(p)、t-GlcA(p)、4,6-GalA(p) | 出现 **Glc**，与组成表的 Man 对不上 |
| 位移归属表（Table 2） | GalA×7、Gal×3、Ara×3、Rha×2、GlcA×1 | 有 Ara / Rha / GlcA，**无 Man** |
| 摘要 | "polymer chain is composed of α-D-galacturonic acid, with branched chains of **α-D-galactose**-α-L-rhamnose, α-D-galactose, and α-L-arabinose" | Table 2 实测为 **β-D-Galp**（H-1 4.39 ppm、C-1 103.29 ppm，与 β 构型一致），摘要写 α **与自表数据矛盾** |

**结论**：该文献的"组成"与"残基/甲基化"是两套口径，且摘要与表数据有冲突。
任何自动化管线都应把这类冲突标 `flagged` 并进入人工复核，而不是静默入库。

---

## 5. 对 ETL 管线的直接改进建议

1. **图片型结构式**（本案例的根本原因）：Fig. 3E 这类"结构以位图呈现"的情形，
   纯文本管线无解。建议在 `table_parser` 增加**图形区域检测 + 结构式 OCR/人工复核队列**：
   先把"疑似结构图"（含 `Galp`/`Rhap`/`→` 字样的图区或引用 Fig. xE 的位置）标记出来，
   产出待复核条目，再决定是否调用 OCR 或人工录入。

2. **补充结论句触发词**（成本最低、收益最高）：在 `extract_linkage_sequence` 之外，
   新增"域级结论"抽取规则，匹配
   `composed of .* (HG|RG-I|RG-II|homogalacturonan|rhamnogalacturonan) domains`、
   `possible structural elements`、`mainly composed of` 等写法，
   产出 **域级结构描述**（写入 `poly_props.backbone/branching` 或新增 `domain_architecture` 字段），
   而不是因为拿不到箭头序列就丢弃整条结论。

3. **新增结构表达等级 `domain_only`**：现有三级
   （`complete` / `repeat_unit` / `composition_only`）无法表达"域级骨架 + 侧链"这一档。
   本案例恰好属于这一档：已知 HG 主链 + RG-I 插入 + 半乳聚糖/阿拉伯聚糖侧链，
   比 `composition_only` 信息量大得多，且**可安全用于域级检索**。

4. **交叉一致性校验**（见评估文档 P2）：组成表 ↔ 残基表 ↔ 甲基化表的组分互查；
   本案例中"Man 只在组成表出现""Ara/Rha/GlcA 只在残基表出现"应自动报警。
