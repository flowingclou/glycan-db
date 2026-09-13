# -*- coding: utf-8 -*-
"""
glycan_etl.py — 糖类数据库数据导入流水线 v3.1（PDF → 结构化入库）

v2 相对 v1 修复四大局限:
  1) 段落分割: 按化合物边界切分文本, 实现"结构↔谱图"逐段正确配对
  2) 2D 谱解析: 新增 HSQC/HMBC/COSY/TOCSY/NOESY/ROESY 相关峰提取
  3) 2D 入库: 相关峰写入 nmr_correlations_2d, 糖苷键证据标记 linkage_evidence
  4) 单糖兜底收紧: 仅段落标题中的糖名作为候选, 不再整篇扫描

v3 新增三大能力（补齐 v2 核心缺口）:
  5) 多糖识别: 重复单元写法 →4)-β-D-Glcp-(1→ / (1→4)-linked β-D-Glcp, sugar_type='poly'
  6) 残基表: 单糖/寡糖/多糖均生成 residues（残基序号/环形式/异头/还原端/连接位点）
  7) 物化+多糖性质: 旋光度/熔点/溶解性/pKa → physicochemical;
      Mn/Mw/PDI/DP/残基摩尔比 → polysaccharide_props

v3.1 修复四项（2026-09-12）:
  8) 独立溶剂/温度提取: "DMSO-d6, 80 °C" 不再整段吞入 solvent,
      solvent='DMSO-d6' + temperature=80.0（1D/2D 头均支持; 无温度则 NULL）
  9) 同核 2D 语义: COSY/TOCSY/NOESY/ROESY 两维均为 ¹H,
      hetero_shift_ppm 置 NULL（不再把 δH 误存为 ¹³C 位移）
  10) 还原端语义: 寡糖还原端为末端残基(如 Glcα1-4Glc 残基2),
      非还原端残基(残基1)标注标题明确构型, 其 C-1 参与糖苷键
  11) R3 误报修复: 寡糖 α/β 混存(不同异头 J)不再强 flag,
      存在与构型参考吻合的异头峰时仅记提示 note

用法:
  python glycan_etl.py --pdf path/to/SI.pdf --dry-run            # 只解析, 不连库
  python glycan_etl.py --pdf path/to/SI.pdf --config cfg.yaml    # 正常入库
  python glycan_etl.py --self-test                                # 内置示例自检

依赖: pip install pdfplumber psycopg2-binary pyyaml
"""
import argparse
import hashlib
import json
import logging
import os
import re
import sys
from dataclasses import dataclass, field, asdict
from typing import List, Optional

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("glycan_etl")

# ============================================================================
# 阶段0: 常量与词典
# ============================================================================

# 常见单糖三字母码（用于结构条目识别）
MONOSACCHARIDES = {
    "Glc", "Gal", "Man", "Fuc", "Xyl", "Rha", "Ara",
    "GlcNAc", "GalNAc", "ManNAc", "Neu5Ac", "Neu5Gc", "Kdo",
    "Fru", "GlcA", "GalA", "IdoA", "MurNAc", "GlcN", "GalN",
}

# 英文全名/俗名 → 三字母码（单糖标题识别）
MONO_SYNONYMS = {
    "Glucose": "Glc", "Galactose": "Gal", "Mannose": "Man",
    "Fucose": "Fuc", "Xylose": "Xyl", "Rhamnose": "Rha",
    "Arabinose": "Ara", "Fructose": "Fru",
    "Glucosamine": "GlcN", "Galactosamine": "GalN",
    "N-Acetylglucosamine": "GlcNAc", "N-acetylglucosamine": "GlcNAc",
    "N-Acetylgalactosamine": "GalNAc", "N-acetylgalactosamine": "GalNAc",
    "N-Acetylneuraminic acid": "Neu5Ac", "N-acetylneuraminic acid": "Neu5Ac",
    "Sialic acid": "Neu5Ac", "sialic acid": "Neu5Ac",
    "Glucuronic acid": "GlcA", "Iduronic acid": "IdoA",
}

# 糖名 → 分子式/分子量（用于规范化, 可扩充）
MONO_FORMULA = {
    "Glc": ("C6H12O6", 180.1559), "Gal": ("C6H12O6", 180.1559),
    "Man": ("C6H12O6", 180.1559), "Fuc": ("C6H12O5", 164.1565),
    "Xyl": ("C5H10O5", 150.1299), "GlcNAc": ("C8H15NO6", 221.2078),
    "GalNAc": ("C8H15NO6", 221.2078), "Neu5Ac": ("C11H19NO9", 309.2699),
    "Fru": ("C6H12O6", 180.1559), "GlcA": ("C6H10O7", 194.1394),
}

# 异头位移合理区间（质控 R1）
ANOMERIC_RANGE = {"13C": (90.0, 110.0), "1H": (4.2, 5.8)}
# J 耦合合理区间（质控 R2, ¹H）
J_RANGE_H1 = (0.5, 18.0)
# 构型-J 参考（质控 R3, 仅异头氢）
ANOMER_J_REF = {"a": 3.8, "b": 7.5}
ANOMER_J_TOL = 2.5

# 正则：IUPAC 连接式 如 Glcα1-4Glc、GalNAcβ1-3Gal
RE_LINKAGE = re.compile(
    r"([A-Z][A-Za-z0-9]+)([αβabAB])(\d+)-(\d+)([A-Z][A-Za-z0-9]+)"
)
# 正则：¹H 峰  δ (多重峰, J=xx Hz, nH, 归属)
RE_PEAK_1H = re.compile(
    r"(?P<shift>\d+\.\d{2,3})\s*\((?P<mult>[sdqtmbr]\w*)\s*,?\s*"
    r"(?:J\s*=\s*(?P<j>[\d.]+)\s*[Hh][Zz])?\s*,?\s*"
    r"(?P<int>\d+\.?\d*)?H?\s*,?\s*"
    r"(?P<assign>[A-Za-z0-9\-', ]*?)\)"
)
# 正则：¹³C 峰  δ (归属)
RE_PEAK_13C = re.compile(
    r"(?P<shift>\d+\.\d{1,2})\s*\(?(?P<assign>[A-Za-z0-9\-', ]*)\)?"
)
# 正则：1D 实验头部 如 "1H NMR (500 MHz, D2O)" 或 "1H NMR (500 MHz, DMSO-d6, 80 °C)"
RE_EXP_HEAD = re.compile(
    r"(?<![A-Za-z])(?P<nucleus>1?[Hh]|1?[Cc])\s*NMR\s*\(?\s*"
    r"(?P<freq>\d+(?:\.\d+)?)\s*[Mm][Hh][Zz]\s*,?\s*"
    r"(?P<solvent>[^(),]+?)(?:\s*,\s*(?P<temp>\d+(?:\.\d+)?)\s*°?\s*[CK])?\s*\)"
)
# 正则：2D 实验头部 如 "HSQC NMR (500 MHz, D2O)"
RE_EXP_2D_HEAD = re.compile(
    r"(?P<exp2d>HSQC|HMBC|COSY|TOCSY|NOESY|ROESY)\s*NMR\s*\(?\s*"
    r"(?P<freq>\d+(?:\.\d+)?)?\s*(?:[Mm]?[Hh][Zz]?)?\s*,?\s*"
    r"(?P<solvent>[^(),]+?)(?:\s*,\s*(?P<temp>\d+(?:\.\d+)?)\s*°?\s*[CK])?\s*\)"
)
# 正则：2D 相关峰  δH/δC (归属) 或  δC/δH (归属)
RE_PEAK_2D = re.compile(
    r"(?P<d1>\d+\.\d{1,3})\s*[/／]\s*(?P<d2>\d+\.\d{1,3})\s*\(?(?P<assign>[^()\n]{0,40})\)?"
)

# ---- 段落边界正则（化合物分块）----
RE_COMPOUND = re.compile(r"(?m)^\s*Compound\s+[\w.\-]+\s*[:：]?")
RE_LINKAGE_HEAD = re.compile(r"(?m)^\s*" + RE_LINKAGE.pattern + r"\s*[:：]?")
RE_MONO_NAME_HEAD = re.compile(
    r"(?mi)^\s*(?:[αβab]-)?[DL]?-?(?:"
    + "|".join(sorted(MONO_SYNONYMS.keys(), key=len, reverse=True))
    + r")\s*[:：]?$"
)
RE_MONO_CODE_HEAD = re.compile(
    r"(?m)^\s*(?:[αβ]-?)?(?:[DL]-?)?(?P<code>"
    + "|".join(sorted(MONOSACCHARIDES, key=len, reverse=True))
    + r")(?:p|f)?(?:[αβ])?\d?\s*[:：]?$"
)

# ============================================================================
# 数据结构
# ============================================================================

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
class Peak2D:
    experiment_2d: str
    proton_shift: float
    hetero_shift: float
    atom_pair: str
    linkage_evidence: bool = False


@dataclass
class NMRExperiment:
    nucleus: str                      # '1H' / '13C' / '2D'
    experiment_2d: Optional[str] = None   # HSQC/HMBC/... (nucleus='2D' 时)
    solvent: str = "CDCl3"
    frequency: Optional[float] = None
    temperature: Optional[float] = None
    ph: Optional[float] = None
    peaks: List[Peak1D] = field(default_factory=list)
    peaks_2d: List[Peak2D] = field(default_factory=list)


@dataclass
class Residue:
    """糖残基（结构组成, 对齐 schema.residues）"""
    residue_seq: int = 1
    monosaccharide_name: str = ""
    ring_form: str = "p"                 # p/f/open/unknown
    anomer: str = "unknown"              # a/b/unknown
    is_reducing_end: bool = False
    parent_carbon: Optional[int] = None  # 糖苷键母体碳位（如 4）
    linkage_branch: int = 0


@dataclass
class Physicochemical:
    """物化性质（对齐 schema.physicochemical）"""
    optical_rotation: Optional[float] = None
    optical_rotation_condition: Optional[str] = None
    melting_point_c: Optional[float] = None
    solubility: Optional[str] = None
    pka: Optional[float] = None


@dataclass
class PolysaccharideProps:
    """多糖专属性质（对齐 schema.polysaccharide_props）"""
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
    """一条待入库的 结构-谱图对（单段落）"""
    sugar_type: str = "mono"
    iupac_short: Optional[str] = None
    glycoct: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    anomer: Optional[str] = None
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


# ============================================================================
# 阶段1: PDF 文本抽取
# ============================================================================

def extract_pdf_text(pdf_path: str) -> str:
    """用 pdfplumber 抽取 PDF 全部文本（按页拼接, 页码可溯源）。"""
    try:
        import pdfplumber
    except ImportError:
        log.error("缺少 pdfplumber, 请先: pip install pdfplumber")
        raise
    pages = []
    with pdfplumber.open(pdf_path) as pdf:
        for i, page in enumerate(pdf.pages, start=1):
            t = page.extract_text() or ""
            pages.append(f"\n<<<PAGE {i}>>>\n{t}")
    return "\n".join(pages)


# ============================================================================
# 阶段2: 结构识别（含单糖标题/同义词）
# ============================================================================

def _anomer_char(ch: str) -> str:
    """把 α/β 或 a/b 统一归一化为 'a'/'b'。"""
    return 'a' if ch in ('α', 'a', 'A') else 'b'


def parse_linkage_iupac(text: str) -> List[str]:
    """从文本中提取 IUPAC 连接式, 如 Glcα1-4Glc → 返回候选结构名。"""
    found = []
    for m in RE_LINKAGE.finditer(text):
        a, anom, c1, c2, b = m.group(1), _anomer_char(m.group(2)), m.group(3), m.group(4), m.group(5)
        if a in MONOSACCHARIDES and b in MONOSACCHARIDES:
            found.append(f"{a}{'α' if anom == 'a' else 'β'}{c1}-{c2}{b}")
    return found


def parse_mono_title(header: str) -> Optional[str]:
    """从段落标题行提取单糖名（三字母码）。

    仅接受"标题行"形式的单糖名：行首、可选 "Compound N:" 前缀、后接
    冒号/行尾/分号。正文中零散出现的糖名（如 "Glc is a common sugar"）
    不会被误识别 —— 这是对"单糖兜底过宽"局限的收紧。
    """
    lines = header.splitlines()[:5]
    for line in lines:
        s = line.strip()
        s = re.sub(r"^Compound\s+[\w.\-]+\s*[:：]\s*", "", s)
        if not s:
            continue
        for syn, code in sorted(MONO_SYNONYMS.items(), key=lambda kv: len(kv[0]), reverse=True):
            if re.match(rf"(?:[αβab]-)?[DL]?-?(?:{syn})(?:p|f)?\s*[:：;,)]?$", s, re.I):
                return code
        for code in sorted(MONOSACCHARIDES, key=len, reverse=True):
            if re.match(rf"(?:[αβ]-?)?(?:[DL]-?)?(?:{code})(?:p|f)?(?:[αβ])?\d?\s*[:：;,)]?$", s):
                return code
    return None


def parse_mono_anomer(header: str) -> Optional[str]:
    """从标题提取异头构型, 如 'β-D-Glcp' → 'b'。"""
    m = re.search(r"([αβaAbB])-", header)
    if m:
        return _anomer_char(m.group(1))
    m = re.search(r"([αβ])(?=\d|$)", header)
    if m:
        return _anomer_char(m.group(1))
    return None


def infer_anomer(iupac: str) -> Optional[str]:
    """从 IUPAC 短名推断异头构型 a/b。"""
    m = re.search(r"[αβ](?=\d)", iupac)
    if m:
        return "a" if m.group(0) == "α" else "b"
    return None


def _merge_formula(f1: str, f2: str) -> str:
    """两个分子式合并为双糖分子式（减去一个 H2O）。如 C6H12O6 + C6H12O6 → C12H22O11。"""
    from collections import defaultdict

    def parse(f):
        d = defaultdict(int)
        for elem, n in re.findall(r"([A-Z][a-z]?)(\d*)", f):
            d[elem] += int(n) if n else 1
        return d

    d = parse(f1)
    for k, v in parse(f2).items():
        d[k] += v
    d["H"] -= 2
    d["O"] -= 1
    return "".join(f"{k}{d[k]}" if d[k] > 1 else k for k in ("C", "H", "O", "N") if d[k] > 0)


def _glycoct_for_mono(code: str, anomer: Optional[str]) -> str:
    """单糖 GlycoCT 示例编码（生产环境建议用 glypy 生成标准编码）。"""
    if anomer:
        return f"RES 1b:{anomer}-{code.lower()}p-1:5|2:x"
    return f"RES 1b:{code.lower()}p-1:5|2:x"


def normalize_sugar_record(iupac: str) -> dict:
    """根据糖名构造 GlycoCT / 分子式 / 构型等。"""
    m = RE_LINKAGE.match(iupac)
    if m:
        a, anom, c1, c2, b = m.group(1), _anomer_char(m.group(2)), m.group(3), m.group(4), m.group(5)
        formula, mw = None, None
        if a in MONO_FORMULA and b in MONO_FORMULA:
            f_a, w_a = MONO_FORMULA[a]
            f_b, w_b = MONO_FORMULA[b]
            formula = _merge_formula(f_a, f_b)
            mw = round(w_a + w_b - 18.015, 4)
        glycoct = f"RES 1b:{anom}-{a.lower()}p-1:5(1:{c2})|2:x,1a:{b.lower()}p-1:5|2:x"
        return {"glycoct": glycoct, "formula": formula, "mw": mw, "anomer": anom}
    # 单糖分支
    code = None
    for c in sorted(MONOSACCHARIDES, key=len, reverse=True):
        if re.search(rf"\b{c}\b", iupac):
            code = c
            break
    if code:
        anom = infer_anomer(iupac) or parse_mono_anomer(iupac)
        formula, mw = MONO_FORMULA.get(code, (None, None))
        return {"glycoct": _glycoct_for_mono(code, anom), "formula": formula, "mw": mw, "anomer": anom}
    return {"glycoct": None, "formula": None, "mw": None, "anomer": None}


# ============================================================================
# 阶段2.5: 多糖重复单元识别 + 残基表生成 + 物化/多糖性质解析（v3 扩展）
# ============================================================================

RE_POLY_UNIT = re.compile(
    r"(?P<parent>\d+)\)\s*-?\s*(?P<anom>[αβabAB])?-?[DL]?-?"
    r"(?P<code>[A-Z][A-Za-z0-9]*?)(?P<ring>[pf])?-\((?P<child>\d+)?[→\-]?"
)
RE_POLY_LINKED = re.compile(
    r"\(1→(?P<parent>\d+)\)\s*-?\s*linked\s+(?P<anom>[αβab])?-?[DL]?-?"
    r"(?P<code>[A-Z][A-Za-z0-9]*?)(?P<ring>[pf])?(?![A-Za-z0-9])",
    re.I,
)
RE_OPT_ROT = re.compile(
    r"\[α\]\s*[Dd]?\s*(?P<temp>\d+(?:\.\d+)?)?\s*[=:]?\s*"
    r"(?P<val>[+-]?\d+(?:\.\d+)?)\s*(?:\((?P<cond>[^)]*)\))?"
)
RE_MP = re.compile(
    r"(?:mp|m\.p\.|melting\s*point)[^0-9]*?(?P<lo>\d+)"
    r"(?:\s*[–\-]\s*(?P<hi>\d+))?\s*°?\s*C", re.I)
RE_PKA = re.compile(r"pKa\s*[=:]?\s*(?P<val>\d+(?:\.\d+)?)", re.I)
RE_MW = re.compile(
    r"(?P<kind>Mn|Mw)\s*[=:]?\s*(?P<val>\d+(?:\.\d+)?)\s*[×xX*]\s*10"
    r"\s*(?:\^)?\s*(?P<exp>\d+)", re.I)
RE_PDI = re.compile(r"(?:PDI|polydispersity)\s*[=:]?\s*(?P<val>\d+(?:\.\d+)?)", re.I)
RE_DP = re.compile(r"(?:DP|degree\s*of\s*polymerization)\s*[=:]?\s*(?P<val>\d+(?:\.\d+)?)", re.I)
RE_RATIO = re.compile(
    r"(?P<a>[A-Z][A-Za-z0-9]*)\s*:\s*(?P<b>[A-Z][A-Za-z0-9]*)\s*[=:]?\s*"
    r"(?P<va>\d+(?:\.\d+)?)\s*:\s*(?P<vb>\d+(?:\.\d+)?)")


def _anomer_greek(a: str) -> str:
    return {"a": "α", "b": "β"}.get(a, "")


def parse_polysaccharide(text: str):
    """识别多糖重复单元写法 → (iupac_short, residues) 或 None。

    支持: →4)-β-D-Glcp-(1→4)-β-D-Glcp-(1→  与  (1→4)-linked β-D-Glcp
    """
    if "→" not in text and "->" not in text:
        return None
    units = [m for m in RE_POLY_UNIT.finditer(text) if m.group("code") in MONOSACCHARIDES]
    if len(units) < 1:
        m = RE_POLY_LINKED.search(text)
        if m and m.group("code") in MONOSACCHARIDES:
            units = [m]
        else:
            return None
    # 连续性检查: 重复单元串需相邻（间隔>40 字符视为分散误配）
    if len(units) > 1:
        for i in range(1, len(units)):
            if units[i].start() - units[i - 1].start() > 40:
                return None
    residues = []
    for i, m in enumerate(units):
        parent = int(m.group("parent")) if m.group("parent") else None
        anom = _anomer_char(m.group("anom")) if m.group("anom") else "unknown"
        ring = m.group("ring") or "p"
        code = m.group("code")
        residues.append(Residue(
            residue_seq=i + 1,
            monosaccharide_name=code,
            ring_form=ring,
            anomer=anom,
            is_reducing_end=False,   # 重复单元无还原端
            parent_carbon=parent,
            linkage_branch=0,
        ))
    # IUPAC 短名: 截取原文重复单元片段（保留前导箭头/中括号）
    start = units[0].start()
    if start > 0 and text[start - 1] in "→[>":
        start -= 1
    iupac = text[start:units[-1].end()]
    iupac = iupac.replace(" ", "").replace("->", "→")
    return iupac, residues


def build_residues(iupac: str, sugar_type: str) -> List[Residue]:
    """单糖/寡糖 → 残基列表（从还原端编号）。"""
    if sugar_type == "mono":
        for c in sorted(MONOSACCHARIDES, key=len, reverse=True):
            if re.search(rf"\b{c}\b", iupac):
                anom = infer_anomer(iupac) or "unknown"
                return [Residue(1, c, ring_form="p", anomer=anom, is_reducing_end=True)]
        return []
    m = RE_LINKAGE.match(iupac)   # oligo: Glcα1-4Glc
    if m:
        a, anom, c1, c2, b = m.group(1), _anomer_char(m.group(2)), m.group(3), m.group(4), m.group(5)
        return [
            # 残基1 = 非还原端（标题明确构型, 其 C-1 参与糖苷键）
            Residue(1, a, ring_form="p", anomer=anom,
                    is_reducing_end=False, parent_carbon=int(c1)),
            # 残基2 = 还原端（α/β 平衡, 构型通常不明确; C-c2 被取代）
            Residue(2, b, ring_form="p", anomer="unknown",
                    is_reducing_end=True, parent_carbon=int(c2)),
        ]
    return []


def _glycoct_for_poly(residues: List[Residue]) -> str:
    """多糖重复单元简化 GlycoCT（生产环境建议 glypy 生成标准编码）。"""
    parts = []
    for r in residues:
        anom = "a" if r.anomer == "a" else "b"
        parts.append(f"{anom}-{r.monosaccharide_name.lower()}{r.ring_form}-1:5"
                     f"(1:{r.parent_carbon or 0})")
    return "B|b:" + ",".join(parts)


def parse_physicochemical(text: str) -> Optional[Physicochemical]:
    """提取旋光度 / 熔点 / 溶解性 / pKa。"""
    pc = Physicochemical()
    m = RE_OPT_ROT.search(text)
    if m:
        pc.optical_rotation = float(m.group("val"))
        cond = (m.group("temp") or "") + (f" {m.group('cond')}" if m.group("cond") else "")
        pc.optical_rotation_condition = cond.strip() or None
    m = RE_MP.search(text)
    if m:
        lo, hi = float(m.group("lo")), m.group("hi")
        pc.melting_point_c = (lo + float(hi)) / 2 if hi else lo
    m = RE_PKA.search(text)
    if m:
        pc.pka = float(m.group("val"))
    m = re.search(r"(?:soluble|solubility)\s*[^.]{0,80}", text, re.I)
    if m:
        pc.solubility = m.group(0).strip()
    if not any([pc.optical_rotation, pc.melting_point_c, pc.solubility, pc.pka]):
        return None
    return pc


def parse_poly_props(text: str) -> Optional[PolysaccharideProps]:
    """提取多糖 Mn/Mw/PDI/DP/残基摩尔比。"""
    pp = PolysaccharideProps()
    for m in RE_MW.finditer(text):
        val = float(m.group("val")) * 10 ** int(m.group("exp"))
        if m.group("kind").upper() == "MN":
            pp.molecular_weight_mn = val
        else:
            pp.molecular_weight_mw = val
    m = RE_PDI.search(text)
    if m:
        pp.polydispersity = float(m.group("val"))
    m = RE_DP.search(text)
    if m:
        pp.degree_of_polymerization = float(m.group("val"))
    m = RE_RATIO.search(text)
    if m and m.group("a") in MONOSACCHARIDES and m.group("b") in MONOSACCHARIDES:
        pp.monosaccharide_ratio = {
            m.group("a"): float(m.group("va")),
            m.group("b"): float(m.group("vb")),
        }
    if not any([pp.molecular_weight_mn, pp.molecular_weight_mw, pp.polydispersity,
                pp.degree_of_polymerization, pp.monosaccharide_ratio]):
        return None
    return pp


# ============================================================================
# 阶段1.5: 段落分割（修复局限1: 结构↔谱图配对）
# ============================================================================

def split_blocks(text: str) -> List[tuple]:
    """按化合物边界把文本切分为 [(header, body), ...]。"""
    markers = []
    for pat in (RE_COMPOUND, RE_LINKAGE_HEAD, RE_MONO_NAME_HEAD, RE_MONO_CODE_HEAD):
        for m in pat.finditer(text):
            markers.append((m.start(), m.group(0).strip()))
    if not markers:
        return [("", text)]
    markers.sort(key=lambda x: x[0])
    # 合并同位置（<30 字符内视为同一标题, 保留更完整的）
    merged = []
    for pos, h in markers:
        if merged and pos - merged[-1][0] < 30:
            # 保留更长的标题
            if len(h) > len(merged[-1][1]):
                merged[-1] = (merged[-1][0], h)
            continue
        merged.append((pos, h))
    blocks = []
    for i, (pos, h) in enumerate(merged):
        end = merged[i + 1][0] if i + 1 < len(merged) else len(text)
        blocks.append((h, text[pos:end]))
    return blocks


# ============================================================================
# 阶段3: NMR 实验解析（1D + 2D, 修复局限2）
# ============================================================================

def _is_linkage_evidence(exp2d: str, assign: str) -> bool:
    """判断 2D 相关峰是否为糖苷键连接证据（HMBC 跨残基相关）。"""
    if exp2d != "HMBC":
        return False
    if "→" in assign or "->" in assign:
        return True
    # 跨残基: 同一相关峰含带撇(如 C-4') 与不带撇(如 H-1) 的编号
    if re.search(r"\d'\s*[,/\-]\s*C-?\d", assign) or re.search(r"C-?\d\s*[,/\-]\s*\d'", assign):
        return True
    # H-1(异头) 与 C-2..6 异位相关
    if re.search(r"[Hh]-?1(?:[αβ]|')?\b", assign) and re.search(r"C-?[2-6](?:'|\b)", assign):
        return True
    return False


def _parse_2d_peaks(seg: str, exp2d: str) -> List[Peak2D]:
    peaks = []
    homonuclear = exp2d in ("COSY", "TOCSY", "NOESY", "ROESY")  # 同核: 两维均为 ¹H
    for pm in RE_PEAK_2D.finditer(seg):
        try:
            d1, d2 = float(pm.group("d1")), float(pm.group("d2"))
        except (TypeError, ValueError):
            continue
        assign = (pm.group("assign") or "").strip()
        if not assign:
            continue
        if homonuclear:
            # 同核 2D: δH-δH, hetero_shift_ppm 语义为 ¹³C, 故留空
            peaks.append(Peak2D(exp2d, d1, None, assign, _is_linkage_evidence(exp2d, assign)))
            continue
        # 顺序判定: 归属中 C 在前则 δC/δH
        m_c = re.search(r"C-?(\d+)", assign)
        m_h = re.search(r"[Hh]-?(\d+)", assign)
        if m_c and m_h and m_c.start() < m_h.start():
            proton, hetero = d2, d1
        else:
            proton, hetero = d1, d2
        peaks.append(Peak2D(exp2d, proton, hetero, assign, _is_linkage_evidence(exp2d, assign)))
    return peaks


def parse_experiments(text: str) -> List[NMRExperiment]:
    """解析 NMR 段落 → 实验对象列表（含 1D 峰与 2D 相关峰）。"""
    heads = [(m.start(), m, "1d") for m in RE_EXP_HEAD.finditer(text)]
    heads += [(m.start(), m, "2d") for m in RE_EXP_2D_HEAD.finditer(text)]
    heads.sort(key=lambda x: x[0])
    exps = []
    for idx, (pos, m, kind) in enumerate(heads):
        end = heads[idx + 1][0] if idx + 1 < len(heads) else pos + 1500
        seg = text[pos:end]
        if kind == "2d":
            exp = NMRExperiment(
                nucleus="2D",
                experiment_2d=m.group("exp2d"),
                solvent=(m.group("solvent") or "").strip() or "CDCl3",
                frequency=float(m.group("freq")) if m.group("freq") else None,
                temperature=float(m.group("temp")) if m.group("temp") else None,
            )
            exp.peaks_2d = _parse_2d_peaks(seg, exp.experiment_2d)
        else:
            nucleus = "13C" if m.group("nucleus").upper() in ("13C", "1C", "C") else "1H"
            exp = NMRExperiment(
                nucleus=nucleus,
                solvent=(m.group("solvent") or "").strip().strip(")") or "CDCl3",
                frequency=float(m.group("freq")) if m.group("freq") else None,
                temperature=float(m.group("temp")) if m.group("temp") else None,
            )
            peak_re = RE_PEAK_1H if nucleus == "1H" else RE_PEAK_13C
            for pm in peak_re.finditer(seg):
                try:
                    shift = float(pm.group("shift"))
                except (TypeError, ValueError):
                    continue
                peak = Peak1D(
                    nucleus=nucleus,
                    shift=shift,
                    multiplicity=pm.group("mult") if nucleus == "1H" else None,
                    j_coupling=float(pm.group("j")) if nucleus == "1H" and pm.group("j") else None,
                    integration=float(pm.group("int")) if nucleus == "1H" and pm.group("int") else None,
                    assignment=(pm.group("assign") or "").strip() or None,
                )
                if peak.assignment and re.search(
                    r"H-?1(?:α|β)?$|C-?1(?:α|β)?$|H-?1[αβ]?[,' ]|C-?1[αβ]?[,' ]|H-1|C-1", peak.assignment
                ):
                    peak.is_anomeric = True
                exp.peaks.append(peak)
        exps.append(exp)
    return exps


# ============================================================================
# 阶段4+5: 配对规范化 + 质控预检
# ============================================================================

def parse_block(header: str, body: str, doi: str, journal: str, year: int) -> List[GlycanRecord]:
    """解析单个化合物段落 → 记录列表（v3: 支持多糖/残基/物化性质）。"""
    text = header + "\n" + body
    sugar_type, iupac, residues, glycoct = "mono", None, [], None
    # 1) 多糖重复单元优先
    poly = parse_polysaccharide(text[:600])
    if poly:
        iupac, residues = poly
        sugar_type = "poly"
        glycoct = _glycoct_for_poly(residues)
    # 2) IUPAC 连接式（寡糖）
    if not iupac:
        cands = parse_linkage_iupac(text[:600])
        if cands:
            iupac, sugar_type = cands[0], "oligo"
    # 3) 单糖标题
    if not iupac:
        iupac = parse_mono_title(header + "\n" + body[:200])
    if not iupac:
        return []
    exps = parse_experiments(text)
    info = normalize_sugar_record(iupac)
    if sugar_type != "poly":
        residues = build_residues(iupac, sugar_type)
    rec = GlycanRecord(
        sugar_type=sugar_type,
        iupac_short=iupac,
        glycoct=glycoct or info["glycoct"],
        molecular_formula=info["formula"],
        molecular_weight=info["mw"],
        anomer=info["anomer"] or parse_mono_anomer(header + "\n" + body[:100]),
        doi=doi, journal=journal, year=year,
        experiments=exps,
        residues=residues,
        physicochemical=parse_physicochemical(text),
        poly_props=parse_poly_props(text) if sugar_type == "poly" else None,
    )
    qc_precheck(rec)
    return [rec]


def qc_precheck(rec: GlycanRecord) -> None:
    """Python 侧质控预检（R1/R2/R3 + 2D 异头碳区间）。"""
    has_2d = any(e.nucleus == "2D" for e in rec.experiments)
    for exp in rec.experiments:
        for pk in exp.peaks:
            rng = ANOMERIC_RANGE.get(pk.nucleus)
            if pk.is_anomeric and rng and not (rng[0] <= pk.shift <= rng[1]):
                rec.qc_notes.append(f"R1: 异头{pk.nucleus}位移 {pk.shift} 越界 {rng}")
                rec.qc_status = "flagged"
            if pk.nucleus == "1H" and pk.j_coupling and not (J_RANGE_H1[0] <= pk.j_coupling <= J_RANGE_H1[1]):
                rec.qc_notes.append(f"R2: J耦合 {pk.j_coupling}Hz 越界")
                rec.qc_status = "flagged"
        for pk in exp.peaks_2d:
            # R1(2D) 仅适用于 HSQC 直接相关, 且 H-1 对应的应是异头碳 C-1
            if pk.experiment_2d == "HSQC" \
                    and re.search(r"[Hh]-?1", pk.atom_pair) \
                    and re.search(r"C-?1", pk.atom_pair) \
                    and not (90.0 <= pk.hetero_shift <= 110.0):
                rec.qc_notes.append(f"R1(2D): HSQC 异头碳 {pk.hetero_shift} 越界区间 90-110")
                rec.qc_status = "flagged"
    # R3: 构型-J 一致性（全局视野: 寡糖常 α/β 混存, 存在吻合峰则不 flag）
    h1_anom_j = [pk.j_coupling for e in rec.experiments for pk in e.peaks
                 if pk.nucleus == "1H" and pk.is_anomeric and pk.j_coupling]
    if h1_anom_j and rec.anomer:
        ref = ANOMER_J_REF.get(rec.anomer)
        if ref:
            mismatched = [j for j in h1_anom_j if abs(j - ref) > ANOMER_J_TOL]
            if mismatched:
                if any(abs(j - ref) <= ANOMER_J_TOL for j in h1_anom_j):
                    rec.qc_notes.append(
                        f"R3(提示): 存在与{rec.anomer}参考 {ref} 不符的异头J {mismatched}（可能α/β混存）")
                else:
                    rec.qc_notes.append(f"R3: {rec.anomer}构型 异头J {mismatched} 均与参考 {ref} 偏差大")
                    rec.qc_status = "flagged"
    if rec.qc_status == "pending":
        rec.qc_status = "passed"
    log.info("QC: %s -> %s (%d notes)", rec.iupac_short, rec.qc_status, len(rec.qc_notes))


# ============================================================================
# 阶段6: 写入 PostgreSQL（1D + 2D, 修复局限3）
# ============================================================================

def _structure_signature(rec: GlycanRecord) -> str:
    """构造可复现的结构指纹（解析不出 GlycoCT 时用作去重键）。"""
    parts = [rec.sugar_type or "?", (rec.iupac_short or "?").strip().lower()]
    if rec.molecular_formula:
        parts.append(rec.molecular_formula.strip())
    if rec.residues:
        parts.append(";".join(
            f"{r.residue_seq}:{r.monosaccharide_name}:{r.ring_form}:{r.anomer}:"
            f"{r.parent_carbon}:{r.linkage_branch}" for r in rec.residues))
    return "|".join(parts)


def resolve_glycoct(rec: GlycanRecord) -> str:
    """返回入库用的结构编码；解析不出 GlycoCT 时生成**确定性占位符**。

    占位符必须逐结构唯一。旧实现回退成空字符串, 于是所有"无 GlycoCT"的
    结构都撞上 sugars.glycoct 的 UNIQUE 约束, 被 ON CONFLICT 互相覆盖:
    多篇文献入库时只有最后一条结构存活, 其余记录的谱图被错误挂到同一条
    结构上(结构↔谱图配对错乱)。占位符含结构指纹, 因此不同结构彼此隔离,
    同一结构重复导入仍能正确命中同一条记录。
    """
    if rec.glycoct and rec.glycoct.strip():
        return rec.glycoct.strip()
    digest = hashlib.sha256(_structure_signature(rec).encode("utf-8")).hexdigest()[:16]
    doi = (rec.doi or "UNKNOWN").strip() or "UNKNOWN"
    return f"UNRESOLVED:{rec.sugar_type or 'unknown'}:{doi}:{digest}"


_COLUMN_CACHE: dict = {}


def _table_columns(cur, table: str) -> set:
    """读取表列名（兼容加过 / 未加 v004 迁移的两种 schema）。"""
    if table not in _COLUMN_CACHE:
        cur.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name = %s",
            (table,),
        )
        _COLUMN_CACHE[table] = {row[0] for row in cur.fetchall()}
    return _COLUMN_CACHE[table]


def insert_record(conn, rec: GlycanRecord) -> None:
    """将一条结构-谱图对写入 v1/v2 表结构（幂等: 重复导入不重复累积）。"""
    cur = conn.cursor()
    cur.execute(
        "INSERT INTO literature (doi, journal, year) VALUES (%s,%s,%s) "
        "ON CONFLICT (doi) DO UPDATE SET journal=EXCLUDED.journal RETURNING source_id",
        (rec.doi, rec.journal, rec.year),
    )
    source_id = cur.fetchone()[0]
    has_2d = any(e.nucleus == "2D" for e in rec.experiments)
    glycoct = resolve_glycoct(rec)
    cur.execute(
        "INSERT INTO sugars (sugar_type, glycoct, glycoct_hash, iupac_short, molecular_formula, "
        "molecular_weight, anomer, structure_confidence, stereochemistry_defined, first_seen_doi) "
        "VALUES (%s,%s, encode(digest(%s,'sha256'),'hex'),%s,%s,%s,%s,%s,%s,%s) "
        # 命中同一结构时只"补空 / 升级", 不覆盖已有的可读结构名
        "ON CONFLICT (glycoct) DO UPDATE SET "
        "  iupac_short = CASE WHEN sugars.iupac_short IS NULL OR sugars.iupac_short = '' "
        "                      THEN EXCLUDED.iupac_short ELSE sugars.iupac_short END, "
        "  structure_confidence = CASE WHEN sugars.structure_confidence = 'confirmed_2d' "
        "                              THEN sugars.structure_confidence "
        "                              ELSE EXCLUDED.structure_confidence END, "
        "  updated_at = now() "
        "RETURNING sugar_id",
        (rec.sugar_type, glycoct, glycoct, rec.iupac_short,
         rec.molecular_formula, rec.molecular_weight, rec.anomer,
         "confirmed_2d" if has_2d else "confirmed_1d",
         rec.anomer is not None, rec.doi),
    )
    sugar_id = cur.fetchone()[0]

    # 幂等: 先清掉本 (结构, 文献) 的旧派生行再重写, 避免重复导入累积脏数据。
    # nmr_shifts_1d / nmr_correlations_2d 由外键 ON DELETE CASCADE 一并清理。
    cur.execute("DELETE FROM nmr_experiments WHERE sugar_id=%s AND source_id=%s",
                (sugar_id, source_id))
    cur.execute("DELETE FROM physicochemical WHERE sugar_id=%s AND source_id=%s",
                (sugar_id, source_id))
    pp_cols = _table_columns(cur, "polysaccharide_props")
    if "source_id" in pp_cols:
        cur.execute("DELETE FROM polysaccharide_props WHERE sugar_id=%s AND source_id=%s",
                    (sugar_id, source_id))
    else:
        cur.execute("DELETE FROM polysaccharide_props WHERE sugar_id=%s", (sugar_id,))

    # v3: residues 残基表（按 (sugar_id, residue_seq) 覆盖）
    for r in rec.residues:
        cur.execute(
            "INSERT INTO residues (sugar_id, residue_seq, monosaccharide_name, ring_form, anomer, "
            "is_reducing_end, parent_carbon, linkage_branch) VALUES (%s,%s,%s,%s,%s,%s,%s,%s) "
            "ON CONFLICT (sugar_id, residue_seq) DO UPDATE SET "
            "  monosaccharide_name=EXCLUDED.monosaccharide_name, ring_form=EXCLUDED.ring_form, "
            "  anomer=EXCLUDED.anomer, is_reducing_end=EXCLUDED.is_reducing_end, "
            "  parent_carbon=EXCLUDED.parent_carbon, linkage_branch=EXCLUDED.linkage_branch",
            (sugar_id, r.residue_seq, r.monosaccharide_name, r.ring_form, r.anomer,
             r.is_reducing_end, r.parent_carbon, r.linkage_branch),
        )
    # v3: physicochemical 物化性质
    if rec.physicochemical:
        pc = rec.physicochemical
        cur.execute(
            "INSERT INTO physicochemical (sugar_id, optical_rotation, optical_rotation_condition, "
            "melting_point_c, solubility, pka, source_id) VALUES (%s,%s,%s,%s,%s,%s,%s)",
            (sugar_id, pc.optical_rotation, pc.optical_rotation_condition,
             pc.melting_point_c, pc.solubility, pc.pka, source_id),
        )
    # v3: polysaccharide_props 多糖专属性质
    if rec.poly_props:
        pp = rec.poly_props
        cols = ["sugar_id", "repeat_unit_formula", "degree_of_polymerization",
                "molecular_weight_mn", "molecular_weight_mw", "polydispersity",
                "monosaccharide_ratio", "backbone", "branching"]
        vals = [sugar_id, pp.repeat_unit_formula, pp.degree_of_polymerization,
                pp.molecular_weight_mn, pp.molecular_weight_mw, pp.polydispersity,
                json.dumps(pp.monosaccharide_ratio) if pp.monosaccharide_ratio else None,
                pp.backbone, pp.branching]
        if "source_id" in pp_cols:
            cols.append("source_id")
            vals.append(source_id)
        placeholders = ",".join(["%s"] * len(vals))
        cur.execute(
            f"INSERT INTO polysaccharide_props ({','.join(cols)}) VALUES ({placeholders})",
            vals,
        )
    for exp in rec.experiments:
        if exp.nucleus == "2D":
            # 2D 实验: 实验头写 nmr_experiments, 相关峰写 nmr_correlations_2d
            cur.execute(
                "INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, "
                "temperature_c, ph, qc_status, assignment_level) "
                "VALUES (%s,%s,'2D',%s,%s,%s,%s,%s,'full') RETURNING experiment_id",
                (sugar_id, source_id, exp.solvent, exp.frequency,
                 exp.temperature, exp.ph, rec.qc_status),
            )
            exp_id = cur.fetchone()[0]
            for pk in exp.peaks_2d:
                cur.execute(
                    "INSERT INTO nmr_correlations_2d (experiment_id, experiment_2d, proton_shift_ppm, "
                    "hetero_shift_ppm, atom_pair, linkage_evidence) "
                    "VALUES (%s,%s,%s,%s,%s,%s)",
                    (exp_id, pk.experiment_2d, pk.proton_shift, pk.hetero_shift,
                     pk.atom_pair, pk.linkage_evidence),
                )
        else:
            cur.execute(
                "INSERT INTO nmr_experiments (sugar_id, source_id, nmr_type, solvent, frequency_mhz, "
                "temperature_c, ph, qc_status, assignment_level) "
                "VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s) RETURNING experiment_id",
                (sugar_id, source_id, exp.nucleus, exp.solvent, exp.frequency,
                 exp.temperature, exp.ph, rec.qc_status,
                 "full" if has_2d else "partial"),
            )
            exp_id = cur.fetchone()[0]
            for pk in exp.peaks:
                cur.execute(
                    "INSERT INTO nmr_shifts_1d (experiment_id, nucleus, shift_ppm, multiplicity, "
                    "j_coupling_hz, integration, assignment_position, is_anomeric) "
                    "VALUES (%s,%s,%s,%s,%s,%s,%s,%s)",
                    (exp_id, pk.nucleus, pk.shift, pk.multiplicity, pk.j_coupling,
                     pk.integration, pk.assignment, pk.is_anomeric),
                )
    conn.commit()
    cur.close()


# ============================================================================
# 主流程
# ============================================================================

def run_etl(pdf_path: str, config: dict, dry_run: bool = False,
            skip_flagged: bool = False):
    log.info("阶段1: 抽取 PDF 文本 %s", pdf_path)
    text = extract_pdf_text(pdf_path)
    doi = config.get("doi", "UNKNOWN")
    journal = config.get("journal")
    year = config.get("year")

    log.info("阶段1.5: 段落分割")
    blocks = split_blocks(text)
    log.info("切分出 %d 个化合物段落", len(blocks))

    records = []
    for header, body in blocks:
        recs = parse_block(header, body, doi, journal, year)
        if recs:
            log.info("段落「%s」→ %s", header[:40] or "(无标题)", recs[0].iupac_short)
            records.extend(recs)
    log.info("共识别 %d 条结构-谱图对", len(records))

    if dry_run:
        out = {"doi": doi, "n_blocks": len(blocks), "records": [asdict(r) for r in records]}
        print(json.dumps(out, ensure_ascii=False, indent=2))
        log.info("DRY-RUN 完成（未连库）")
        return records

    import psycopg2
    conn = psycopg2.connect(**config["db"])
    n, skipped = 0, 0
    for rec in records:
        if rec.qc_status == "flagged" and skip_flagged:
            log.warning("跳过违规记录 %s: %s", rec.iupac_short, rec.qc_notes)
            skipped += 1
            continue
        if rec.qc_status == "flagged":
            # 默认保留可疑数据入库并打上 flagged 标记: 对真值库而言,
            # 丢弃比保留风险更高(无法回溯), 由下游查询自行决定是否采信。
            log.warning("入库但标记 flagged: %s: %s", rec.iupac_short, rec.qc_notes)
        insert_record(conn, rec)
        n += 1
    conn.close()
    log.info("入库完成, 共写入 %d 条（跳过 %d 条）", n, skipped)


# ============================================================================
# 自检
# ============================================================================

SAMPLE_TEXT = """
Compound 1: Glcα1-4Glc
1H NMR (500 MHz, D2O) δ 5.41 (d, J = 3.8 Hz, 1H, H-1'), 5.23 (d, J = 3.8 Hz, 1H, H-1),
4.66 (d, J = 7.9 Hz, 1H, H-1β), 3.85 (m, 2H, H-6), 3.52 (dd, J = 9.9, 3.8 Hz, 1H, H-2).
13C NMR (126 MHz, D2O) δ 100.6 (C-1'), 96.8 (C-1), 78.3 (C-4), 73.9 (C-3), 72.5 (C-2).
HSQC NMR (500 MHz, D2O) δ 5.41/100.6 (H-1', C-1'), 3.52/72.5 (H-2, C-2).
HMBC NMR (500 MHz, D2O) δ 5.41/78.3 (H-1', C-4).
Compound 2: β-D-Glcp
1H NMR (500 MHz, D2O) δ 4.64 (d, J = 7.9 Hz, 1H, H-1), 3.52 (dd, J = 9.9, 7.9 Hz, 1H, H-2).
13C NMR (126 MHz, D2O) δ 96.8 (C-1).
"""


def self_test():
    text = SAMPLE_TEXT
    log.info("SELF-TEST: 段落分割")
    blocks = split_blocks(text)
    assert len(blocks) == 2, f"期望 2 段落, 实际 {len(blocks)}"
    h1, h2 = blocks[0][0], blocks[1][0]
    log.info("OK: 段落1=「%s」 段落2=「%s」", h1, h2)

    log.info("SELF-TEST: 段落1 结构↔谱图配对")
    recs1 = parse_block(*blocks[0], "10.xx", "J. Org. Chem.", 2020)
    assert len(recs1) == 1
    r1 = recs1[0]
    assert r1.iupac_short == "Glcα1-4Glc", r1.iupac_short
    assert r1.sugar_type == "oligo"
    exps = r1.experiments
    assert len(exps) == 4, f"期望 4 实验(1H/13C/HSQC/HMBC), 实际 {len(exps)}"
    hsqc = next(e for e in exps if e.experiment_2d == "HSQC")
    hmbc = next(e for e in exps if e.experiment_2d == "HMBC")
    assert len(hsqc.peaks_2d) == 2, f"HSQC 期望 2 峰, 实际 {len(hsqc.peaks_2d)}"
    assert len(hmbc.peaks_2d) == 1
    assert hmbc.peaks_2d[0].linkage_evidence is True, "HMBC H-1'→C-4 应为糖苷键证据"
    log.info("OK: 4 实验, HSQC %d 峰, HMBC %d 峰, 糖苷键证据=%s",
             len(hsqc.peaks_2d), len(hmbc.peaks_2d), hmbc.peaks_2d[0].linkage_evidence)

    log.info("SELF-TEST: 段落2 单糖识别 + 构型")
    recs2 = parse_block(*blocks[1], "10.xx", "J. Org. Chem.", 2020)
    assert len(recs2) == 1
    r2 = recs2[0]
    assert r2.iupac_short == "Glc", r2.iupac_short
    assert r2.anomer == "b", r2.anomer
    assert r2.sugar_type == "mono"
    log.info("OK: 单糖 %s (anomer=%s)", r2.iupac_short, r2.anomer)

    log.info("SELF-TEST: 单糖兜底收紧（正文零散糖名不产生记录）")
    # 只有 NMR 峰文本、无结构标题的输入 → 0 记录
    stray = "Glc is a common sugar. 1H NMR (500 MHz, D2O) d 5.22 (d, J=3.8 Hz, 1H, H-1)."
    recs3 = parse_block("", stray, "10.xx", "J", 2020)
    assert len(recs3) == 0, f"正文零散糖名不应产生记录, 实际 {len(recs3)}"
    log.info("OK: 无结构标题 → 0 记录")

    log.info("SELF-TEST: 多糖重复单元识别 + 残基表")
    poly_text = ("→4)-β-D-Glcp-(1→4)-β-D-Glcp-(1→\n"
                 "1H NMR (500 MHz, D2O) δ 4.52 (br s, 1H, H-1), 3.31 (m, 1H, H-2).\n"
                 "13C NMR (126 MHz, D2O) δ 103.8 (C-1), 79.5 (C-4).")
    rp = parse_block("", poly_text, "10.xx", "Carbohydr. Polym.", 2023)
    assert len(rp) == 1, f"多糖应识别, 实际 {len(rp)}"
    assert rp[0].sugar_type == "poly", rp[0].sugar_type
    assert len(rp[0].residues) == 2, f"重复单元应含 2 残基, 实际 {len(rp[0].residues)}"
    assert rp[0].residues[0].parent_carbon == 4, rp[0].residues[0].parent_carbon
    assert rp[0].residues[0].anomer == "b"
    log.info("OK: 多糖 %s, %d 残基, 首残基 anomer=%s parent_carbon=%s",
             rp[0].iupac_short, len(rp[0].residues), rp[0].residues[0].anomer,
             rp[0].residues[0].parent_carbon)

    log.info("SELF-TEST: (1→4)-linked 单残基多糖写法")
    poly_linked = ("Cellulose, (1→4)-linked β-D-Glcp\n"
                   "1H NMR (500 MHz, D2O) δ 4.52 (br s, 1H, H-1).")
    rl = parse_block("", poly_linked, "10.xx", "Biomacromolecules", 2023)
    assert len(rl) == 1, f"(1→4)-linked 应识别, 实际 {len(rl)}"
    assert rl[0].sugar_type == "poly" and rl[0].residues[0].monosaccharide_name == "Glc"
    assert rl[0].residues[0].parent_carbon == 4
    log.info("OK: %s → %d 残基 parent_carbon=%s",
             rl[0].iupac_short, len(rl[0].residues), rl[0].residues[0].parent_carbon)

    log.info("SELF-TEST: 寡糖残基表")
    r_oligo = build_residues("Glcα1-4Glc", "oligo")
    assert len(r_oligo) == 2, f"寡糖应含 2 残基, 实际 {len(r_oligo)}"
    assert not r_oligo[0].is_reducing_end and r_oligo[0].anomer == "a", "残基1应为非还原端 α"
    assert r_oligo[1].is_reducing_end and r_oligo[1].parent_carbon == 4, "残基2应为还原端, C-4 被取代"
    log.info("OK: Glcα1-4Glc → 2 残基, 残基1(α,非还原端) 残基2(还原端,C4被取代)")

    log.info("SELF-TEST: 物化性质 + 多糖性质")
    phys_text = ("[α]D20 +50.2 (c 1.0, H2O); mp 210-212 °C; pKa 3.5. "
                 "Mw 1.2 × 10^5, Mn 8.5 × 10^4, PDI 1.4, DP 20, Glc:Gal = 7:3.")
    pc = parse_physicochemical(phys_text)
    assert pc and pc.optical_rotation == 50.2 and abs(pc.melting_point_c - 211.0) < 0.01 and pc.pka == 3.5
    pp = parse_poly_props(phys_text)
    assert pp and pp.molecular_weight_mw == 120000 and pp.molecular_weight_mn == 85000
    assert pp.polydispersity == 1.4 and pp.degree_of_polymerization == 20
    assert pp.monosaccharide_ratio == {"Glc": 7.0, "Gal": 3.0}
    log.info("OK: 旋光度 %.1f, 熔点 %.1f°C, pKa %.1f, Mw %g, Mn %g, PDI %.2f, DP %.0f, 比例 %s",
             pc.optical_rotation, pc.melting_point_c, pc.pka, pp.molecular_weight_mw,
             pp.molecular_weight_mn, pp.polydispersity, pp.degree_of_polymerization,
             pp.monosaccharide_ratio)

    log.info("SELF-TEST: 全部通过 ✓")


def main():
    ap = argparse.ArgumentParser(description="糖类数据库 PDF→库 ETL v2")
    ap.add_argument("--pdf", help="SI 补充材料 PDF 路径")
    ap.add_argument("--config", default="etl_config.yaml", help="配置文件")
    ap.add_argument("--dry-run", action="store_true", help="只解析不连库")
    ap.add_argument("--self-test", action="store_true", help="内置示例自检")
    ap.add_argument("--skip-flagged", action="store_true",
                    help="质控违规记录不入库（默认入库并标记 flagged，便于回溯）")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return
    if not args.pdf:
        ap.error("需要 --pdf 或 --self-test")

    config = {}
    if os.path.exists(args.config):
        import yaml
        with open(args.config, encoding="utf-8") as f:
            config = yaml.safe_load(f) or {}
    run_etl(args.pdf, config, dry_run=args.dry_run,
            skip_flagged=args.skip_flagged)


if __name__ == "__main__":
    main()
