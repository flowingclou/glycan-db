#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
glycan_table_parser.py -- 多糖结构表征文献「表格式数据」解析器
================================================================
定位：补齐 glycan_etl_v3.py 对「以表格形式呈现的多糖结构表征文献」
识别为 0 条的缺口。glycan_etl_v3 面向段落式（正文叙述）数据；本解析
器面向「分子量/组成表 + ¹H/¹³C 化学位移归属表 + 正文键连描述」三类
表格型内容，输出与 glycan_etl_v3 风格一致的结构化 dry-run 结果
（不连数据库，控制台 JSON 输出）。

解析能力（三类）：
  1) 分子量与物化组成表（Mp/Mw/Mn/RT + 各单糖 mol%）
  2) ¹H/¹³C 化学位移归属表（行列式：残基编码 × H-1~H-6a/6b /
     C-1~C-6 / -OMe，H 行与 C 行成对出现）
  3) 正文中甲基化/糖苷键连描述（如 "mainly 4-GalA(p)"、
     "→4)-α-D-GalpA-(1→"）

设计要点：
  - 复用 glycan_etl_v3.py 的 GlycanRecord / Residue / NMRExperiment /
    Peak1D / PolysaccharideProps 数据类与 RE_POLY_UNIT、MONOSACCHARIDES
    等正则/词典，保证条目输出风格统一。
  - 表格以「词级坐标聚类」重建（pdfplumber extract_words），解决
    双栏 PDF 提取文本时表格与正文交错、以及无框线对齐式表格
    无法被 extract_tables() 识别的问题。
  - 化学位移表按「表头列锚点最近邻」把每个数值归位到 *exact*
    位置（H-1..H-6a/6b / C-1..C-6 / -OMe），不依赖文本顺序。

用法：
  python3 glycan_table_parser.py <polysaccharide_paper.pdf> [--json out.json]
输出：
  控制台打印统计 + dry-run JSON；可 --json 落盘。
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import OrderedDict, defaultdict
from dataclasses import asdict
from typing import Any, Dict, List, Optional, Tuple

# 将仓库根加入 sys.path，使包内相对引用可用（支持以脚本方式直接运行）
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# ----------------------------------------------------------------------------
# 复用 glycan_etl.core（原 glycan_etl_v3.py）：数据类 + 正则/词典（保持项目风格统一）
# ----------------------------------------------------------------------------
try:
    from glycan_etl import core as etl

    GlycanRecord = etl.GlycanRecord
    Residue = etl.Residue
    Peak1D = etl.Peak1D
    NMRExperiment = etl.NMRExperiment
    Physicochemical = etl.Physicochemical
    PolysaccharideProps = etl.PolysaccharideProps
    MONOSACCHARIDES = set(etl.MONOSACCHARIDES)
    RE_POLY_UNIT = etl.RE_POLY_UNIT
    ETL_VER = getattr(etl, "__version__", "v3")
    _ANOM_K = etl._anomer_greek
except Exception as exc:  # pragma: no cover - 兜底，避免 import 失败即崩溃
    print(f"[warn] 未能导入 glycan_etl_v3，使用内置简化数据类: {exc}", file=sys.stderr)
    from dataclasses import dataclass, field

    @dataclass
    class Peak1D:
        nucleus: str
        shift: float
        multiplicity: Optional[str] = None
        j_coupling: Optional[float] = None
        integration: Optional[float] = None
        assignment: Optional[str] = None
        is_anomeric: bool = False

    @dataclass
    class NMRExperiment:
        nucleus: str
        solvent: str = "D2O"
        frequency: Optional[float] = None
        temperature: Optional[float] = None
        ph: Optional[float] = None
        peaks: List[Peak1D] = field(default_factory=list)

    @dataclass
    class Residue:
        residue_code: str = ""
        residue_seq: int = 1
        monosaccharide_name: str = ""
        ring_form: str = "p"
        anomer: str = "unknown"
        is_reducing_end: bool = False
        parent_carbon: Optional[int] = None
        linkage_branch: int = 0

    @dataclass
    class Physicochemical:
        optical_rotation: Optional[float] = None
        optical_rotation_condition: Optional[str] = None
        melting_point_c: Optional[float] = None
        solubility: Optional[str] = None
        pka: Optional[float] = None

    @dataclass
    class PolysaccharideProps:
        repeat_unit_formula: Optional[str] = None
        degree_of_polymerization: Optional[float] = None
        molecular_weight_mn: Optional[float] = None
        molecular_weight_mw: Optional[float] = None
        polydispersity: Optional[float] = None
        monosaccharide_ratio: Optional[dict] = None
        backbone: Optional[str] = None
        branching: Optional[str] = None

    @dataclass
    class GlycanRecord:
        sugar_type: str = "poly"
        iupac_short: Optional[str] = None
        mol_name: Optional[str] = None
        doi: Optional[str] = None
        journal: Optional[str] = None
        year: Optional[int] = None
        nmr_page: Optional[int] = None
        residues: List[Residue] = field(default_factory=list)
        physicochemical: Optional[Physicochemical] = None
        poly_props: Optional[PolysaccharideProps] = None
        experiments: List[NMRExperiment] = field(default_factory=list)
        qc_notes: List[str] = field(default_factory=list)
        qc_status: str = "pending"

    @dataclass
    class _ETL:
        MONOSACCHARIDES = {
            "Glc", "Gal", "Man", "Fuc", "Xyl", "Rha", "Ara",
            "GlcNAc", "GalNAc", "ManNAc", "Neu5Ac", "Neu5Gc", "Kdo",
            "Fru", "GlcA", "GalA", "IdoA", "MurNAc", "GlcN", "GalN",
        }

    MONOSACCHARIDES = set(_ETL.MONOSACCHARIDES)
    RE_POLY_UNIT = None
    ETL_VER = "standalone-fallback"
    _ANOM_K = lambda a: {"a": "α", "b": "β"}.get(a, "")

PARSER_VER = "1.0.0"

# ----------------------------------------------------------------------------
# 基础工具：词级坐标 → 行重建
# ----------------------------------------------------------------------------
ROW_DY = 5          # 同一视觉行的 top 容差（px）
ND_TOKENS = {"–", "-", "—", "nd", "n.d.", "n.d", "…"}


def pdf_pages(pdf_path: str) -> List[Dict[str, Any]]:
    """逐页提取：文本 + 词级坐标；并将词聚类为带坐标的『行』。"""
    try:
        import pdfplumber
    except ImportError as e:
        raise SystemExit(
            f"[fatal] 需要 pdfplumber（pip install pdfplumber）: {e}")

    pages: List[Dict[str, Any]] = []
    with pdfplumber.open(pdf_path) as pdf:
        n_pages = len(pdf.pages)
        for i, page in enumerate(pdf.pages, 1):
            words = page.extract_words(
                keep_blank_chars=False, use_text_flow=False,
                extra_attrs=["size", "fontname"])
            lines = cluster_rows(words, dy=ROW_DY)
            text = page.extract_text() or ""
            pages.append({
                "page": i, "n_pages": n_pages, "text": text,
                "words": words, "lines": lines,
                "width": page.width, "height": page.height,
            })
    return pages


def cluster_rows(words: List[dict], dy: int = ROW_DY) -> List[dict]:
    """把单词按 top 聚成『行』；行内按 x0 排序。返回 [{top, words, text}]。"""
    rows: List[dict] = []
    for w in sorted(words, key=lambda w: (w["top"], w["x0"])):
        best, bd = None, dy + 1
        for r in rows:
            d = abs(r["top"] - w["top"])
            if d < bd:
                best, bd = r, d
        if best is None or best["top"] == 0 and rows:  # pragma: no cover
            pass
        if best is None:
            best = {"top": w["top"], "words": []}
            rows.append(best)
        best["words"].append(w)
        # 更新行为该簇的『顶部』（取最小值），保持后续归簇稳定
        best["top"] = min(best["top"], w["top"])
    for r in rows:
        r["words"] = sorted(r["words"], key=lambda w: (w["x0"]))
        r["text"] = " ".join(w["text"] for w in r["words"])
    rows.sort(key=lambda r: r["top"])
    return rows


def _is_num(s: str) -> bool:
    s = s.replace(",", "")
    try:
        float(s)
        return True
    except ValueError:
        pass
    return bool(re.match(r"^\s*-?\d+(?:\.\d+)?\s*/\s*-?\d+(?:\.\d+)?\s*$", s))


def _num(s: str) -> Optional[float]:
    s = s.replace(",", "").strip()
    try:
        return float(s)
    except ValueError:
        pass
    # 'a/b' 并存信号（如 H-5 差向 5.01/5.06）取均值
    m = re.match(r"^\s*(-?\d+(?:\.\d+)?)\s*/\s*(-?\d+(?:\.\d+)?)\s*$", s)
    if m:
        return (float(m.group(1)) + float(m.group(2))) / 2.0
    # 仅取首个数字
    m = re.match(r"^(-?\d+(?:\.\d+)?)", s)
    return float(m.group(1)) if m else None


# ----------------------------------------------------------------------------
# 1) 分子量与物化组成表（键-值两列表）
# ----------------------------------------------------------------------------
# 标签 → 规范化字段（单位括号如 (Da)/(min)/(%) 可选）
RE_KV_LABEL = re.compile(
    r"^\s*(?P<label>"
    r"(?P<rt>RT|retention\s*time)\s*\(min\)"
    r"|(?P<mp>Mp|M\(p\)|peak\s+molecular\s+weight|molecular\s+weight\s*\(Mp\))\s*\([^)]*\)?"
    r"|(?P<mw>Mw|M\(w\)|weight-average\s+molecular\s+weight)\s*\([^)]*\)?"
    r"|(?P<mn>Mn|M\(n\)|number-average\s+molecular\s+weight)\s*\([^)]*\)?"
    r"|(?P<gravity>GI|AI|AA|IC50|EC50)\s*\)?"
    r"|(?P<sugar>[A-Z][A-Za-z]{2,}[A-Za-z0-9]*)\s*\(%\)"
    r")\s*$"
)
RE_CAPTION = re.compile(r"^Table\s+(\d+)\b", re.I)


def parse_molecular_kv_table(page: Dict[str, Any],
                             start_idx: int,
                             max_rows: int = 60) -> Tuple[dict, int]:
    """从第 start_idx 个『行』开始解析键值两列表（分子量/组成）。

    返回 (data_dict, next_index)。
    数据行形如：`RT (min)  32.955`、`Mw (Da)  168,471`、`GalA (%)  85.10`。
    与正文交错的干扰列（双栏排版右栏）通过 x0 窗口过滤。
    """
    lines = page["lines"]
    kv: "OrderedDict[str, str]" = OrderedDict()
    last_table = start_idx
    skipped = 0
    for idx in range(start_idx, min(start_idx + max_rows, len(lines))):
        ln = lines[idx]
        # 当前行中所有词（整页），按 x 位置取左侧表区（<0.5*page_w）
        ws = sorted(ln["words"], key=lambda w: w["x0"])
        label_words = [w for w in ws if w["x0"] < 160]
        val_words = [w for w in ws if 160 <= w["x0"] <= 330]
        label_txt = " ".join(w["text"] for w in label_words).strip()
        val_txt = " ".join(w["text"] for w in val_words).strip()

        m = RE_KV_LABEL.match(label_txt)
        if m and val_txt:
            # 值只取第一个数值 token（避免同列混入额外文本）
            val = " ".join(t for t in val_txt.split() if _is_num(t))
            if val:
                key = (m.group("rt") and "RT" or
                       m.group("mp") and "Mp" or
                       m.group("mw") and "Mw" or
                       m.group("mn") and "Mn" or
                       m.group("sugar"))
                kv[key] = val.split()[0]
                last_table = idx
                skipped = 0
                continue
        # 表头行（Sample / CPP 之类）只作提示，不算数据；不重置跳过计数
        if re.search(r"Sample", label_txt) or re.match(r"^\s*CPP\s*$", label_txt):
            skipped = 0
            continue
        # 连续空/非表行则退出
        if not label_txt and not val_txt:
            skipped += 1
            if skipped >= 3:
                break
            last_table = idx
            continue
        # 出现下一张表/图的题注则退出
        if RE_CAPTION.match(label_txt) or re.match(r"^Fig(ure)?\.?\s*\d", label_txt):
            break
        skipped += 1
        if skipped >= 4:
            break
        last_table = idx
    return kv, last_table + 1


def normalize_molecular_kv(kv: Dict[str, str]) -> Dict[str, Any]:
    """把原始键-值对规范化为结构化字段。"""
    out: Dict[str, Any] = {}
    ratio: Dict[str, float] = {}
    for k, v in kv.items():
        f = _num(v)
        if k in ("RT", "Mp", "Mw", "Mn"):
            out[k if k == "RT" else f"molecular_weight_{k.lower()}"] = f
        else:  # 单糖 mol%
            ratio[k] = f
    if ratio:
        out["monosaccharide_mol_pct"] = ratio
    return out


# ----------------------------------------------------------------------------
# 2) ¹H/¹³C 化学位移归属表（行列式：残基 × 位置）
# ----------------------------------------------------------------------------
# 位置标题行元素：`1 2 3 4 5 6a/6b -OMe`
RE_POS_HEADER = re.compile(r"^(?=.*\b1\b)(?=.*\b2\b)(?=.*\b3\b)(?=.*\b4\b)"
                           r"(?=.*\b5\b)(?=.*(?:6a/6b|-OMe|OMe))",
                           re.I | re.X)


def _extract_residue_code(tokens: List[dict]) -> str:
    """Residues 列（x<90）的文本拼接；如 'GE1,4' / '4MeGlcA' / 'GAt'。"""
    ws = [w for w in tokens if w["x0"] < 90]
    return "".join(w["text"] for w in ws).strip()


def _extract_iupac_desc(tokens: List[dict]) -> str:
    """Residues 列与 H/C 列之间（88<=x<195）的 IUPAC 描述。"""
    ws = [w for w in tokens if 88 <= w["x0"] < 195]
    return " ".join(w["text"] for w in ws).strip()


def parse_shift_table(page: Dict[str, Any], start_idx: int,
                      max_rows: int = 80) -> Tuple[List[dict], int]:
    """解析化学位移归属表（依赖词坐标列对齐 + H/C 成对行）。

    返回 (residue_entries, next_index)。每个 entry:
      {code, iupac, parent_carbon, anomer, reducing_end,
       positions: {pos: {"H": float|None, "C": float|None}}, rows}
    """
    lines = page["lines"]
    entries: List[dict] = []
    # ---- 定位位置标题行，得到列锚点（pos -> x_center） ----
    col_anchor: "OrderedDict[str, float]" = OrderedDict()
    header_idx = None
    for idx in range(start_idx, min(start_idx + 6, len(lines))):
        ws = sorted(lines[idx]["words"], key=lambda w: w["x0"])
        txts = [w["text"].strip() for w in ws if w["x0"] > 190]
        joined = " ".join(txts)
        if re.search(r"6a/6b", joined) or re.search(r"-OMe", joined) \
                and re.search(r"\b5\b", joined):
            pos_names = [t for t in txts if t in
                         ("1", "2", "3", "4", "5", "6", "6a/6b", "-OMe", "OMe")]
            if len(pos_names) >= 4:
                nxt = 0
                for w in ws:
                    if w["text"] in pos_names:
                        col_anchor[w["text"]] = (w["x0"] + w["x1"]) / 2.0
                header_idx = idx
                break
    if not col_anchor:
        return [], start_idx

    # fallback：仅数字标题（无 6a/6b）/ -OMe 也行
    if len(col_anchor) < 4:
        return [], start_idx

    # H/C 标记列锚点（'H' / 'C' 所在 x 的近似值，用于区分 H 行与 C 行）
    hc_anchor = None
    for idx in range(header_idx, min(header_idx + 3, len(lines))):
        for w in lines[idx]["words"]:
            if w["text"] in ("H", "C") and 190 <= w["x0"] <= 215:
                hc_anchor = (w["x0"] + w["x1"]) / 2.0
                break
        if hc_anchor:
            break
    if hc_anchor is None:
        hc_anchor = 205.0

    anchor_list = list(col_anchor.items())

    def col_of(x_center: float) -> Optional[str]:
        best, bd = None, 1e9
        for name, ax in anchor_list:
            d = abs(x_center - ax)
            if d < bd:
                best, bd = name, d
        return best if bd <= 30 else "?"

    # ---- 遍历行，识别 残基H行 / C行 ----
    h_row_meta: list = []   # (idx, code, iupac, {pos:val})
    next_idx = start_idx
    rows_consumed = 0
    i = start_idx
    while i < len(lines) and i < start_idx + max_rows:
        ln = lines[i]
        ws = sorted(ln["words"], key=lambda w: w["x0"])
        code = _extract_residue_code(ws)
        desc = _extract_iupac_desc(ws)
        # 判定 H/C 标记
        mark = None
        for w in ws:
            if abs((w["x0"] + w["x1"]) / 2.0 - hc_anchor) <= 12 and w["text"] in ("H", "C"):
                mark = w["text"]
                break
        nums = [w for w in ws if _is_num(w["text"]) and w["x0"] >= 220]
        is_data = bool(code) and mark == "H" and len(nums) >= 2

        if is_data:
            pos: "OrderedDict[str, dict]" = OrderedDict(
                (n, {"H": None, "C": None}) for n in col_anchor)
            for w in nums:
                col = col_of((w["x0"] + w["x1"]) / 2.0)
                if col and col in pos and pos[col]["H"] is None:
                    pos[col]["H"] = _num(w["text"])
            h_row_meta.append({"idx": i, "code": code, "iupac": desc,
                               "pos": pos})
            rows_consumed = 0
            i += 1
            continue

        # C 行：紧随上一残基 H 行（top 相近），有 'C' 标记 + 数值
        if mark == "C" and h_row_meta and not code and len(nums) >= 1:
            prev = h_row_meta[-1]
            if abs(ln["top"] - lines[prev["idx"]]["top"]) <= 25:
                for w in nums:
                    col = col_of((w["x0"] + w["x1"]) / 2.0)
                    if col and col in prev["pos"] and prev["pos"][col]["C"] is None:
                        prev["pos"][col]["C"] = _num(w["text"])
                rows_consumed = 0
                i += 1
                continue

        # 描述行（含占位符 – ）挂到上一残基
        if h_row_meta and desc and not code and mark is None:
            prev = h_row_meta[-1]
            if abs(ln["top"] - lines[prev["idx"]]["top"]) <= 15:
                # 仅当已记录 iupac 为空时补充
                if not prev["iupac"]:
                    prev["iupac"] = desc
                rows_consumed = 0
                i += 1
                continue

        rows_consumed += 1
        if rows_consumed >= 4:
            break
        i += 1

    next_idx = i

    # ---- 组装 entries ----
    for hr in h_row_meta:
        pos_out: "OrderedDict[str, dict]" = OrderedDict()
        for name, p in hr["pos"].items():
            if p["H"] is not None or p["C"] is not None:
                pos_out[name] = {"H": p["H"], "C": p["C"]}
        entries.append({
            "code": hr["code"], "iupac": hr["iupac"], "positions": pos_out,
            "row_index": hr["idx"],
        })
    return entries, next_idx


# 从 IUPAC 描述提取单糖块：`α -D-Gal p A -6-OMe-(1→` → seg=Gal, ring=p, acid=A
RE_MONO_SEG = re.compile(
    r"(?:[αβab])?\s*-?\s*[DL]?-?\s*"
    r"(?P<seg>[A-Z][A-Za-z0-9]*?)\s*(?P<ring>[pf])?\s*(?P<acid>A)?"
    r"\s*(?:-|\(|$)"
)
RE_IUPAC_POS = re.compile(r"(?:→\s*|\(1→\s*)(?P<pos>\d+(?:\s*,\s*\d+)*)\)")
RE_OMe_SEG = re.compile(r"(?:-?\d+\s*-?\s*O\s*-?\s*Me|[-–]\s*(?:\d+-)?OMe|OMe)")


def parse_mono_from_iupac(iupac: str, code: str = "",
                          MONOS: Optional[set] = None) -> dict:
    """从 IUPAC 描述/残基编码解析单糖身份。

    返回 {"mono": str, "ring": str, "anomer": str, "parent": int|None,
           "branch": list[int], "oMe": int|None, "terminal": bool}
    """
    MONOS = MONOS or MONOSACCHARIDES
    out = {"mono": "", "ring": "p", "anomer": "unknown",
           "parent": None, "branch": [], "oMe": None, "terminal": False}
    if not (iupac or code):
        return out

    m = re.search(r"([αβ])", iupac)
    if m:
        out["anomer"] = "a" if m.group(1) == "α" else "b"

    mp = RE_IUPAC_POS.search(iupac)
    if mp:
        poss = [int(x) for x in re.findall(r"\d+", mp.group("pos"))]
        out["parent"] = poss[0]
        if len(poss) > 1:
            out["branch"] = poss[1:]
    else:
        # 无 →pos：视作端基 t-
        out["terminal"] = True

    msub = RE_OMe_SEG.search(iupac)
    if msub:
        nm = re.search(r"(\d+)\s*-?\s*(?:O\s*-?\s*Me|OMe)", iupac)
        out["oMe"] = int(nm.group(1)) if nm else None

    # 糖名：剥离 O-Me 取代段后匹配 IUPAC 主体
    sugar_src = RE_OMe_SEG.sub("", iupac)
    ms = RE_MONO_SEG.search(sugar_src)
    if ms:
        seg, ring, acid = ms.group("seg"), ms.group("ring"), ms.group("acid")
        cand = seg + (acid or "")
        if cand in MONOS:
            out["mono"] = cand
        elif seg in MONOS:
            out["mono"] = seg
        out["ring"] = ring or "p"
        if out["mono"]:
            return out
    # 退化 IUPAC → 修正段后再试（如 "p A" 拆词粘连）
    for k in sorted(MONOS, key=len, reverse=True):
        if k in iupac or (k[-1] in "A" and k[:-1] in iupac and "-p" in iupac):
            out["mono"] = k
            break
    if out["mono"]:
        return out
    # 最终退化：从编码前缀推断（GE1,4→GalA；4MeGlcA→GlcA；Rα→Rha）
    CODE_HINT = {"GE": "GalA", "GA": "GalA", "GAt": "GalA", "T": "GalA",
                 "G": "Gal", "R": "Rha", "RA": "Rha", "At": "Ara",
                 "A": "Ara", "Rha": "Rha", "Ara": "Ara", "Gal": "Gal",
                 "Glc": "Glc", "Man": "Man"}
    code_clean = re.sub(r"[αβab,]", "", code)
    m2 = re.match(r"(?:(\d+)-?Me)?[_-]?([A-Za-z][A-Za-z0-9]*?)(?:,\d+)*$",
                  code_clean)
    if m2:
        head = m2.group(2)
        # 直接命中已知简写
        for pre, full in CODE_HINT.items():
            if head == pre or head.startswith(pre):
                out["mono"] = full
                break
        if not out["mono"]:
            for cand in (head + "A", head, head.rstrip("A")):
                if cand in MONOS:
                    out["mono"] = cand
                    break
        if not out["mono"] and m2.group(1):
            # 如 '4MeGlcA' 的 Me 段已被吞入 group(2)? 异常处理
            pass
        if not out["mono"]:
            head_l = head.lower()
            for k in sorted(MONOS, key=len, reverse=True):
                if k.lower().startswith(head_l) or head_l.startswith(k.lower()):
                    out["mono"] = k
                    break
    return out


def interpret_residue_entry(entry: dict, seq: int) -> Residue:
    """把归属表条目映射为 Residue（含糖名/环/异头/母体碳位推断）。"""
    code = entry["code"]
    iupac = entry["iupac"] or ""
    pk = parse_mono_from_iupac(iupac, code)

    r = Residue()
    r.residue_seq = seq
    r.monosaccharide_name = pk["mono"] if pk["mono"] else code.split(",")[0]
    r.ring_form = pk["ring"]
    r.anomer = pk["anomer"]
    r.is_reducing_end = bool(re.search(r"R[αβab]$", code)) or \
        "reducing" in iupac.lower()
    r.parent_carbon = pk["parent"]
    if pk["branch"]:
        r.linkage_branch = pk["branch"][0]
    return r


# ----------------------------------------------------------------------------
# 3) 正文键连描述（甲基化糖 + IUPAC 重复单元）
# ----------------------------------------------------------------------------
# 甲基化残基：t-GalA(p)、4-Glc(p)、2,4-GalA(p)、4,6-GalA(p)、3,6-β-D-Galp
RE_METH_R = re.compile(
    r"(?P<pos>\d+(?:\s*[,，]\s*\d+)*|t|T)\s*-\s*"
    r"(?P<sugar>[A-Z][A-Za-z0-9]*)\s*\((?P<ring>[pf])\)",
    re.I,
)
# IUPAC 重复单元（复用 etl.RE_POLY_UNIT 思路，含 t- 端基）
RE_UNIT_R = re.compile(
    r"(?:→|\(1→)\s*(?P<pos>\d+(?:,\d+)*)\)\s*-?\s*(?:OMe-)?"
    r"(?P<anom>[αβab])?-?[DL]?-?"
    r"(?P<code>[A-Z][A-Za-z0-9]*?)(?P<ring>[pf])?-?\("
)


def entry_to_repeat_linkages(entry: dict) -> List[dict]:
    """从化学位移归属表的残基条目生成 repeat_unit 键连（统一走单体解析）。

    条目 iupac 形如：
      →4)-α-D-GalpA-6-OMe-(1→        → 4-GalA(p,a) [+ 6-OMe 取代]
      →3,4)-α-D-GalpA-(1→            → 3,4-GalA(p,a)
      →4)-β-D-GalpA                  → 4-Gal(p,b)   (还原端 Rβ)
      α-L-Araf-(1→  /  α-D-GalpA-(1→ → t-Ara(f,a) / t-GalA(p,a) (端基)
    """
    code = entry.get("code") or ""
    iupac = entry.get("iupac") or code
    pk = parse_mono_from_iupac(iupac, code)
    links: List[dict] = []

    if not pk["mono"]:
        return links
    anomer = pk["anomer"]
    if pk["terminal"]:
        links.append({"linkage_type": "repeat_unit", "position": "t",
                      "sugar": pk["mono"], "ring": pk["ring"],
                      "anomer": anomer, "count": 1, "src": code})
    else:
        pos = [str(pk["parent"])] + [str(b) for b in pk["branch"]]
        links.append({"linkage_type": "repeat_unit",
                      "position": ",".join(pos),
                      "sugar": pk["mono"], "ring": pk["ring"],
                      "anomer": anomer, "count": 1, "src": code})
    if pk["oMe"] is not None:
        links.append({"linkage_type": "repeat_unit_substituent",
                      "position": str(pk["oMe"]), "sugar": "OMe",
                      "ring": "", "anomer": "unknown", "count": 1, "src": code})
    return links


def _iupac_anomer(s: str) -> str:
    m = re.search(r"([αβ])", s)
    return "a" if m and m.group(1) == "α" else ("b" if m and m.group(1) == "β"
                                                else "unknown")


def parse_linkages(full_text: str, residue_entries: Optional[List[dict]] = None
                   ) -> List[dict]:
    """从全文提取键连列表：
       'methylation' —— 甲基化分析给出的取代位置（t-GalA(p) 等）
       'repeat_unit' —— IUPAC 呈现键连（→4)-α-D-GalpA-(1→ 等）
    去重（保留首次出现）。
    """
    links: "OrderedDict[str, dict]" = OrderedDict()

    def _add(l: dict):
        k = (f"{l['linkage_type']}:{l.get('position','?')}-{l.get('sugar','?')}"
             f"-{l.get('ring','')}-{l.get('anomer','')}")
        if k in links:
            links[k]["count"] += 1
            return
        links[k] = l

    # IUPAC 重复单元：优先来自归属表条目（信息最完整、最准）
    if residue_entries:
        for ent in residue_entries:
            for l in entry_to_repeat_linkages(ent):
                _add(l)
    else:
        for m in RE_UNIT_R.finditer(full_text):
            sugar = m.group("code")
            if sugar not in MONOSACCHARIDES:
                continue
            _add({"linkage_type": "repeat_unit",
                  "position": re.sub(r"\s+", "", m.group("pos")),
                  "sugar": sugar,
                  "ring": (m.group("ring") or "p"),
                  "anomer": ("a" if m.group("anom") == "α" else
                             "b" if m.group("anom") == "β" else "unknown"),
                  "count": 1})

    # 甲基化残基：全文扫描
    for m in RE_METH_R.finditer(full_text):
        sugar = m.group("sugar")
        if sugar not in MONOSACCHARIDES:
            continue
        _add({"linkage_type": "methylation",
              "position": re.sub(r"\s+", "", m.group("pos")),
              "sugar": sugar,
              "ring": (m.group("ring") or "p"),
              "anomer": "unknown",
              "count": 1})
    return list(links.values())


# ----------------------------------------------------------------------------
# Dry-run 组装
# ----------------------------------------------------------------------------
def _poly_summary_name(residues: List[Residue]) -> str:
    """文献未给出可读多糖名时，用残基组成构造可区分的结构标识。

    例: 16 残基的山楂多糖 → ``poly[GalA x6,Gal x3,Ara x3,Rha x2,GlcA]``。
    这比统一占位串 "unknown polysaccharide" 有用得多：下游 AI 平台回查 /
    展示时至少能区分不同结构，并与 core.resolve_glycoct() 的结构指纹互证。
    """
    counts: "OrderedDict[str, int]" = OrderedDict()
    for r in residues:
        name = (r.monosaccharide_name or "").strip()
        if name:
            counts[name] = counts.get(name, 0) + 1
    if not counts:
        return "unresolved polysaccharide"
    parts = ",".join(
        f"{k} x{v}" if v > 1 else k
        for k, v in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))
    return f"poly[{parts}]"


def build_record(meta: dict, kv_norm: dict, entries: List[dict],
                 linkages: List[dict], nmr_page: int) -> GlycanRecord:
    rec = GlycanRecord()
    rec.sugar_type = "poly"
    rec.iupac_short = meta.get("polysaccharide") or kv_norm.get("name")
    rec.mol_name = meta.get("polysaccharide")
    rec.doi = meta.get("doi")
    rec.journal = meta.get("journal")
    rec.year = meta.get("year")
    rec.nmr_page = nmr_page

    props = PolysaccharideProps()
    if "molecular_weight_mw" in kv_norm:
        props.molecular_weight_mw = kv_norm["molecular_weight_mw"]
    if "molecular_weight_mn" in kv_norm:
        props.molecular_weight_mn = kv_norm["molecular_weight_mn"]
    if "monosaccharide_mol_pct" in kv_norm:
        ratio = kv_norm["monosaccharide_mol_pct"]
        props.monosaccharide_ratio = {k: v for k, v in ratio.items()}

    # backbone / branching 由键连推导（简单拼接）；Mp/RT 由 molecular_table 段承载
    meth = [l for l in linkages if l["linkage_type"] == "methylation"]
    if meth:
        props.branching = ", ".join(
            f"{l['position']}-{l['sugar']}({l['ring']})" for l in meth)
    unit = [l for l in linkages if l["linkage_type"] == "repeat_unit"]
    if unit:
        props.backbone = "; ".join(
            f"{l.get('position','?')}-{l['sugar']}({l['ring']},{l['anomer']})"
            for l in unit)
    rec.poly_props = props

    # 残基 + 实验（两实验：1H / 13C）
    h_exp = NMRExperiment(nucleus="1H", solvent="D2O")
    c_exp = NMRExperiment(nucleus="13C", solvent="D2O")
    for seq, ent in enumerate(entries, 1):
        res = interpret_residue_entry(ent, seq)
        rec.residues.append(res)
        code = ent["code"]
        for posname, pv in ent["positions"].items():
            if pv["H"] is not None:
                anom = "H-1" == _label(posname, "H")
                h_exp.peaks.append(Peak1D(
                    nucleus="1H", shift=pv["H"], assignment=f"{code} {posname}",
                    is_anomeric=anom))
            if pv["C"] is not None:
                anom = "C-1" == _label(posname, "C")
                c_exp.peaks.append(Peak1D(
                    nucleus="13C", shift=pv["C"], assignment=f"{code} {posname}",
                    is_anomeric=anom))
    if h_exp.peaks:
        rec.experiments.append(h_exp)
    if c_exp.peaks:
        rec.experiments.append(c_exp)

    # 文献没给出可读多糖名时，用残基组成构造可区分的结构标识，
    # 避免所有多糖都写成同一个 "unknown polysaccharide" 占位串。
    if not rec.iupac_short or str(rec.iupac_short).strip().lower() in (
            "unknown polysaccharide", "unknown", "unresolved polysaccharide"):
        rec.iupac_short = _poly_summary_name(rec.residues)

    rec.qc_notes.append("table-parser: molecular table (Table-like 1)")
    rec.qc_notes.append(f"table-parser: shift table with {len(entries)} residues")
    rec.qc_status = "dry-run-ok"
    return rec


def _label(pos: str, kind: str) -> str:
    if pos in ("1", "2", "3", "4", "5"):
        return f"{kind}-{pos}"
    if "6" in pos:
        return f"{kind}-6a/6b"
    if "OMe" in pos or "Me" in pos:
        return f"{kind}-OCH3"
    return f"{kind}-{pos}"


# ----------------------------------------------------------------------------
# 元数据（标题/期刊/年份/DOI/多糖名）
# ----------------------------------------------------------------------------
JOURNAL_PAT = re.compile(
    r"(International\s+Journal\s+of\s+Biological\s+Macromolecules"
    r"|Carbohydrate\s+Polymers|Food\s+Hydrocolloids"
    r"|Journal\s+of\s+Agricultural\s+and\s+Food\s+Chemistry"
    r"|Food\s+Chemistry|Carbohydrate\s+Research)",
    re.I)

# 多糖名提取要排除的"非名称"词（避免 'of a polysaccharide' → 'of'）
_NAME_STOPWORDS = {
    "a", "an", "the", "of", "and", "or", "in", "on", "for", "from", "this",
    "that", "its", "new", "novel", "crude", "one", "two", "various", "several",
    "structural", "characterisation", "characterization", "structure", "study",
    "extraction", "purification", "analysis", "effects", "effect", "role",
}


def _despace(s: str) -> str:
    """还原被字间距拆散的文本（Elsevier 页眉常抽成 'J o u r n a l'）。"""
    return re.sub(r"(?<=\b[A-Za-z])\s(?=[A-Za-z]\b)", "", s)


def extract_title(pages: List[dict]) -> Optional[str]:
    """按「最大字号 + 位于页面上半部」取标题行。

    旧实现取「首页第一条长度 > 40 的文本行」：在 Elsevier 版式上会稳定
    命中页眉的期刊名（期刊名被字间距拆开后恰好超长），于是一篇多糖论文
    的 title 变成 "International Journal of Biological Macromolecules"，
    并连带使多糖名提取（依赖 title 中的 'X polysaccharide'）一起失效。
    """
    if not pages:
        return None
    pg = pages[0]
    words = pg.get("words") or []
    if not words:
        return None
    sizes = [w.get("size") or 0 for w in words]
    mx = max(sizes) if sizes else 0
    if mx <= 0:
        return None
    height = pg.get("height") or max((w.get("bottom") or 0) for w in words)
    head = [w for w in words
            if (w.get("size") or 0) >= mx - 0.6
            and (w.get("top") or 0) < height * 0.55]
    if not head:
        return None
    rows = cluster_rows(head, dy=ROW_DY)
    if not rows:
        return None
    cand = max(rows, key=lambda r: len(r["text"]))
    txt = re.sub(r"\s+", " ", cand["text"]).strip()
    if len(txt) < 25 or JOURNAL_PAT.search(txt) or re.search(
            r"contents|elsevier|springer|wiley|volume\s+\d+|https?://|@", txt, re.I):
        return None
    return txt or None


def extract_meta(pages: List[dict]) -> dict:
    first = pages[0]["text"] if pages else ""
    meta = {"doi": None, "journal": None, "year": None, "title": None,
            "polysaccharide": None}
    m = re.search(r"10\.\d{4,9}/[^\s,;\]]+", first)
    if m:
        meta["doi"] = m.group(0).rstrip(".")
    # 期刊名可能被字间距拆散，先还原再匹配
    m = JOURNAL_PAT.search(_despace(first)) or JOURNAL_PAT.search(first)
    if m:
        meta["journal"] = re.sub(r"\s+", " ", m.group(1)).strip()
    m = re.search(r"\(?(20\d{2})\)?", first)
    if m:
        meta["year"] = int(m.group(1))
    meta["title"] = extract_title(pages)
    if not meta["title"]:
        # 兜底：旧的长行启发式
        lines = [l.strip() for l in first.splitlines() if l.strip()]
        for l in lines[1:]:
            if len(l) > 40 and not re.search(r"@|href|Contents", l):
                meta["title"] = l
                break
    # 多糖名：标题中紧邻 'polysaccharide' 的那个词（如 'hawthorn polysaccharide'）
    if meta["title"]:
        m = re.search(r"([A-Za-z][\w-]*)\s+polysaccharide", meta["title"], re.I)
        if m:
            cand = m.group(1).strip()
            if cand.lower() not in _NAME_STOPWORDS:
                meta["polysaccharide"] = cand
    return meta


# ----------------------------------------------------------------------------
# 主流程
# ----------------------------------------------------------------------------
def run(pdf_path: str, verbose: bool = True) -> dict:
    if not os.path.isfile(pdf_path):
        raise SystemExit(f"[fatal] 文件不存在: {pdf_path}")
    print(f"[info] 解析 PDF: {pdf_path}")
    pages = pdf_pages(pdf_path)
    meta = extract_meta(pages)
    print(f"       页数={pages[0]['n_pages']}  "
          f"doi={meta['doi']}  journal={meta['journal']}  year={meta['year']}")

    kv_all: "OrderedDict[str, str]" = OrderedDict()
    entries_all: List[dict] = []
    nmr_page = None
    kv_page = None
    kvcaption: List[str] = []
    shiftcaption: List[str] = []

    for page in pages:
        ln = page["lines"]
        for idx, line in enumerate(ln):
            txt = line["text"].strip()
            mcap = RE_CAPTION.match(txt)
            if not mcap:
                continue
            # 有题注：往后探测是『分子量表』还是『化学位移表』
            cap_no = mcap.group(1)
            probe = " ".join(l["text"] for l in ln[idx:idx + 8])
            if re.search(r"weight|composition|physicochem|physical", probe, re.I) \
                    and not re.search(r"chemical shift|attribution|residu", probe, re.I):
                kv, ni = parse_molecular_kv_table(page, idx + 1)
                if kv:
                    kvcaption.append(f"Table {cap_no}")
                    kv_page = page["page"]
                    kv_all.update(kv)
            elif re.search(r"chemical shift|attribution|C\s*chemical|residu",
                           probe, re.I):
                ents, ni = parse_shift_table(page, idx + 1, max_rows=90)
                if ents:
                    shiftcaption.append(f"Table {cap_no}")
                    nmr_page = page["page"]
                    entries_all.extend(ents)

    # 正文全文（用于键连）：按页序拼接自然文本
    full_text = "\n".join(
        f"\n===== PAGE {pg['page']} =====\n{pg['text']}" for pg in pages)
    linkages = parse_linkages(full_text, residue_entries=entries_all)

    kv_norm = normalize_molecular_kv(kv_all)
    kv_norm["name"] = meta.get("polysaccharide") or "unknown polysaccharide"
    rec = build_record(meta, kv_norm, entries_all, linkages, nmr_page)

    # 计算位移条数
    n_shifts = sum(len(e.peaks) for e in rec.experiments)
    stats = {
        "pages": pages[0]["n_pages"],
        "n_polysaccharide_records": 1 if rec.residues or kv_norm else 0,
        "n_residues": len(rec.residues),
        "n_shifts": n_shifts,
        "n_linkages": len(linkages),
        "n_kv_rows": len(kv_all),
        "molecular_table": kv_page,
        "shift_table": nmr_page,
        "kv_captions": kvcaption,
        "shift_captions": shiftcaption,
    }

    d = {
        "parser": f"glycan_table_parser v{PARSER_VER}",
        "etl_backend": ETL_VER,
        "document": {
            "file": os.path.abspath(pdf_path),
            "filename": os.path.basename(pdf_path),
            "pages": pages[0]["n_pages"],
        },
        "meta": meta,
        "stats": stats,
        "molecular_table": kv_norm,
        "records": [asdict(rec)],
        "linkages": linkages,
    }
    if verbose:
        print("\n========== 解析统计 ==========")
        for k, v in stats.items():
            print(f"  {k}: {v}")
        print(f"\n[info] 键连明细 ({len(linkages)}):")
        for l in linkages:
            print(f"  - [{l['linkage_type']}] "
                  f"{l['position']}-{l['sugar']}({l['ring']}) "
                  f"x{l['count']} anomer={l['anomer']}")
    return d


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="glycan_table_parser")
    ap.add_argument("pdf", help="多糖结构表征文献 PDF 路径")
    ap.add_argument("--json", help="可选：输出 JSON 到文件")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    d = run(args.pdf, verbose=not args.quiet)
    blob = json.dumps(d, ensure_ascii=False, indent=2)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            f.write(blob)
        print(f"\n[info] JSON 已写入: {os.path.abspath(args.json)}")
    else:
        print("\n========== DRY-RUN JSON ==========")
        print(blob)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception as e:
        print(f"[fatal] {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(1)
