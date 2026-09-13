# -*- coding: utf-8 -*-
"""交叉一致性校验层（cross-consistency checks）。

背景
----
同一篇文献里，"单糖组成表"、"位移归属表（残基）"、"甲基化分析"三者应当是
互相印证的。实测却在真实文献中频繁出现互不相容：

* 组成百分比合计远小于 100%（表被截断 / 漏抽行）；
* 某单糖只出现在组成表（如 Man），却不在任何残基表里，反之亦然；
* 甲基化分析出现残基表没有的糖（如 Glc）；
* 摘要写的构型（α-D-galactose）与自家归属表（β-D-Galp）矛盾。

这些冲突此前**完全静默**：数据照常入库，下游拿到的是互相打架的"真值"。
本模块把它们变成显式的 QC 发现（finding），随记录一起落库/进报告。

规则编号沿用仓库的 R 系列，避免与 v003 的 R1-R10 冲突，这里用 C 系列：

  C1 组成百分比合计异常（应为 ~100%）
  C2 单糖种类在"组成表 ↔ 残基表"之间不匹配
  C3 单糖种类在"残基表 ↔ 甲基化分析"之间不匹配
  C4 分子式 ↔ 残基组成式不匹配（仅结构完整时做严格比对）
  C5 摘要/正文构型 与 归属表构型冲突
  C6 甲基化给出的连接类型 与 位移表残基类型 无法对应（数量级差异）
"""
from __future__ import annotations

import re
from collections import Counter, defaultdict
from typing import Any, Dict, Iterable, List, Optional

# 单糖别名 → 标准码（把 "Galacturonic acid"/"GalA" 等统一）
_MONO_ALIASES = {
    "galacturonicacid": "GalA", "galacturonic": "GalA", "gala": "GalA",
    "ge": "GalA", "ga": "GalA", "gat": "GalA",
    "glucuronicacid": "GlcA", "glucuronic": "GlcA", "glca": "GlcA",
    "galactose": "Gal", "glucose": "Glc", "mannose": "Man",
    "rhamnose": "Rha", "arabinose": "Ara", "fucose": "Fuc",
    "xylose": "Xyl", "ribose": "Rib", "fructose": "Fru",
    "glucosamine": "GlcN", "galactosamine": "GalN",
    "nacetylglucosamine": "GlcNAc", "nacetylgalactosamine": "GalNAc",
    "glucuronicacidresidue": "GlcA",
}

# 分子式里每种单糖残基的原子数（游离糖形式；成苷时每键脱一分子水）
_MONO_ATOMS = {
    "Glc": {"C": 6, "H": 12, "O": 6}, "Gal": {"C": 6, "H": 12, "O": 6},
    "Man": {"C": 6, "H": 12, "O": 6}, "All": {"C": 6, "H": 12, "O": 6},
    "Alt": {"C": 6, "H": 12, "O": 6}, "Gul": {"C": 6, "H": 12, "O": 6},
    "Ido": {"C": 6, "H": 12, "O": 6}, "Tal": {"C": 6, "H": 12, "O": 6},
    "Fuc": {"C": 6, "H": 12, "O": 5}, "Rha": {"C": 6, "H": 12, "O": 5},
    "Qui": {"C": 6, "H": 12, "O": 5},
    "Xyl": {"C": 5, "H": 10, "O": 5}, "Ara": {"C": 5, "H": 10, "O": 5},
    "Rib": {"C": 5, "H": 10, "O": 5}, "Lyx": {"C": 5, "H": 10, "O": 5},
    "GlcA": {"C": 6, "H": 10, "O": 7}, "GalA": {"C": 6, "H": 10, "O": 7},
    "ManA": {"C": 6, "H": 10, "O": 7}, "IdoA": {"C": 6, "H": 10, "O": 7},
    "GulA": {"C": 6, "H": 10, "O": 7},
    "GlcNAc": {"C": 8, "H": 15, "N": 1, "O": 6},
    "GalNAc": {"C": 8, "H": 15, "N": 1, "O": 6},
    "ManNAc": {"C": 8, "H": 15, "N": 1, "O": 6},
    "GlcN": {"C": 6, "H": 13, "N": 1, "O": 5},
    "GalN": {"C": 6, "H": 13, "N": 1, "O": 5},
    "Fru": {"C": 6, "H": 12, "O": 6},
}


def canon_mono(name: str) -> Optional[str]:
    """把各种写法（GalA / Galacturonic acid / gal a / GE1,4）归一到标准单糖码。

    返回 None 表示无法识别为单糖。
    """
    if not name:
        return None
    s = str(name).strip()
    # 残基编码（如 GE1,4 / 4MeGlcA / Rha1,2）先剥掉位置与取代前缀
    s = re.sub(r"^\d+(?:,\d+)*-?\s*", "", s)
    s = re.sub(r"^\d+Me", "", s)
    # 异头构型前缀：希腊字母可省连字符；ASCII a/b 必须带连字符，
    # 否则会把糖名尾字母吃掉（GalA→GlA、Gal→Gl），导致全部识别失败。
    s = re.sub(r"^(?:[αβ]\s*-?\s*[DL]?\s*-?|[abAB]\s*-\s*[DL]?\s*-?)", "", s)
    s = re.sub(r"[\s\-_,]", "", s)
    s = re.sub(r"[pf](?=[A-Z]|$)", "", s)          # Galp → Gal
    s = re.sub(r"\d+(?:,\d+)*$", "", s)            # Rha1,2 → Rha
    if not s:
        return None
    low = s.lower()
    if low in _MONO_ALIASES:
        return _MONO_ALIASES[low]
    for code in sorted(_MONO_ATOMS, key=len, reverse=True):
        if low == code.lower():
            return code
    # 退化：包含匹配。必须按候选长度降序，否则 GlcNAc 会被 GlcA 抢先命中。
    for alias, code in sorted(_MONO_ALIASES.items(), key=lambda kv: -len(kv[0])):
        if alias in low:
            return code
    for code in sorted(_MONO_ATOMS, key=len, reverse=True):
        if code.lower() in low:
            return code
    return None


def parse_formula(formula: str) -> Optional[Dict[str, int]]:
    """解析分子式字符串（如 ``C12H22O11``）→ 元素计数字典。"""
    if not formula:
        return None
    m = re.findall(r"([A-Z][a-z]?)(\d*)", str(formula))
    if not m:
        return None
    out: Dict[str, int] = defaultdict(int)
    for elem, n in m:
        if not elem:
            continue
        out[elem] += int(n) if n else 1
    return dict(out) if out else None


def _attr(obj: Any, key: str, default=None):
    """统一读取 dataclass 属性或 dry-run dict 字段。"""
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def composition_from_residues(residues: Iterable[Any]) -> Counter:
    """残基列表 → 单糖计数（用于与组成式/分子式对账）。

    同时支持 ``core.Residue`` 与 dry-run 出来的 dict（键名相同）。
    """
    c: Counter = Counter()
    for r in residues or []:
        code = canon_mono(_attr(r, "monosaccharide_name", "") or "")
        if code:
            c[code] += 1
    return c


def formula_from_composition(composition: Iterable[Any],
                             n_bonds: Optional[int] = None) -> Optional[Dict[str, int]]:
    """由残基计数推算分子式。

    ``n_bonds`` 缺省取 ``残基数 - 1``（一条链）；每形成一个糖苷键脱一分子水。
    """
    # 组成式字符串（如 "GalA7,Ara3,Glc2"）必须先按逗号拆项；否则会被当成
    # 字符序列逐字符遍历，永远解析不出任何残基（实测踩过）。
    if isinstance(composition, str):
        items: List[Any] = [x for x in re.split(r"[,\s]+", composition) if x]
    else:
        items = list(composition or [])

    counts: Counter = Counter()
    total_units = 0
    for item in items:
        if isinstance(item, str):
            m = re.fullmatch(r"\s*([A-Za-z]+)\s*(\d*)\s*", item)
            if not m:
                continue
            code = canon_mono(m.group(1))
            n = int(m.group(2)) if m.group(2) else 1
        else:                                    # (name, count) 元组
            code = canon_mono(str(item[0]))
            n = int(item[1])
        if not code or code not in _MONO_ATOMS:
            continue
        counts[code] += n
        total_units += n
    if not counts:
        return None
    atoms: Dict[str, int] = defaultdict(int)
    for code, n in counts.items():
        for elem, cnt in _MONO_ATOMS[code].items():
            atoms[elem] += cnt * n
    bonds = (total_units - 1) if n_bonds is None else int(n_bonds)
    atoms["H"] -= 2 * bonds
    atoms["O"] -= bonds
    return {k: v for k, v in atoms.items() if v}


def _fmt_atoms(atoms: Dict[str, int]) -> str:
    order = ["C", "H", "N", "O"]
    parts = [f"{e}{atoms[e]}" if atoms.get(e, 0) > 1 else e
             for e in order if atoms.get(e, 0)]
    for e in sorted(set(atoms) - set(order)):
        parts.append(f"{e}{atoms[e]}" if atoms[e] > 1 else e)
    return "".join(parts)


def check_composition_percent_sum(percentages: Dict[str, Any],
                                  tol: float = 5.0) -> List[Dict[str, Any]]:
    """C1：组成百分比合计应接近 100%。"""
    if not percentages:
        return []
    try:
        total = sum(float(v) for v in percentages.values() if v is not None)
    except (TypeError, ValueError):
        return []
    if abs(total - 100.0) > tol:
        return [{
            "rule": "C1", "severity": "flagged",
            "message": (f"单糖组成百分比合计 {total:.2f}%，偏离 100% 超过 {tol:g}%"
                        f"（疑似表格漏抽行或口径不一致）"),
            "fields": {"metric": "mol_pct_sum", "total": round(total, 3),
                       "expected": 100.0, "tolerance": tol},
        }]
    return []


def check_component_sets(percentages: Dict[str, Any],
                         residues: Iterable[Any],
                         methylation: Optional[Iterable[Any]] = None
                         ) -> List[Dict[str, Any]]:
    """C2/C3：组成表 ↔ 残基表 ↔ 甲基化分析 的单糖种类是否对得上。"""
    findings: List[Dict[str, Any]] = []
    pct_codes = {canon_mono(k) for k in (percentages or {})}
    pct_codes.discard(None)
    res_codes = {c for c in composition_from_residues(residues)}
    if pct_codes and res_codes:
        only_pct = sorted(pct_codes - res_codes)
        only_res = sorted(res_codes - pct_codes)
        if only_pct or only_res:
            findings.append({
                "rule": "C2", "severity": "flagged",
                "message": ("单糖种类在'组成表 ↔ 残基表'之间不匹配："
                            f"仅组成表有 {only_pct or '—'}；仅残基表有 {only_res or '—'}"),
                "fields": {"metric": "component_set",
                           "only_in_percentages": only_pct,
                           "only_in_residues": only_res},
            })
    if methylation:
        meth_codes = {canon_mono(x) for x in methylation}
        meth_codes.discard(None)
        if meth_codes and res_codes:
            only_m = sorted(meth_codes - res_codes)
            if only_m:
                findings.append({
                    "rule": "C3", "severity": "flagged",
                    "message": (f"甲基化分析出现位移归属表未列出的单糖 {only_m}"
                                "（两套数据口径不一致，需人工确认）"),
                    "fields": {"metric": "methylation_vs_residues",
                               "only_in_methylation": only_m},
                })
    return findings


def check_formula_vs_residues(formula: Optional[str],
                              composition: Optional[Any],
                              structure_level: Optional[str] = None,
                              n_bonds: Optional[int] = None,
                              tol: float = 1.0) -> List[Dict[str, Any]]:
    """C4：文献给出的分子式 ↔ 由残基组成推算的分子式。

    只在能明确对应时做严格比对：结构不完整（composition_only / domain_only）
    时，聚合度未知，换算出的分子式没有可比性，仅在有明确残基计数时提示。
    """
    lit = parse_formula(formula) if formula else None
    if not lit or not composition:
        return []
    calc = formula_from_composition(composition, n_bonds=n_bonds)
    if not calc:
        return []
    diff = {e: calc.get(e, 0) - lit.get(e, 0)
            for e in set(calc) | set(lit) if calc.get(e, 0) != lit.get(e, 0)}
    level_note = ""
    if structure_level in ("composition_only", "domain_only"):
        level_note = "（结构不完整，聚合度未知，差异可能来自聚合度而非解析错误）"
        severity = "warning"
    else:
        severity = "flagged" if diff else "info"
    if not diff:
        return []
    return [{
        "rule": "C4", "severity": severity,
        "message": (f"分子式对账不符：文献 {formula}，按残基组成推算 "
                    f"{_fmt_atoms(calc)}，差异 {diff}{level_note}"),
        "fields": {"metric": "formula_vs_residues",
                   "literature": formula, "calculated": _fmt_atoms(calc),
                   "diff": {k: v for k, v in diff.items()}},
    }]


RE_ABSTRACT_ANOMER = re.compile(
    r"\b(?P<anom>[αβ])\s*-\s*[DL]?\s*-?\s*(?P<sugar>galactose|glucose|mannose|"
    r"rhamnose|arabinose|fucose|xylose)\b", re.I)


def check_anomer_conflict(text: str, residues: Iterable[Any]) -> List[Dict[str, Any]]:
    """C5：正文/摘要陈述的构型 与 归属表实测构型冲突。

    只在**同一种单糖**上比较（如摘要说 α-D-galactose，而表里该糖残基全是
    β-D-Galp），避免把不同残基的构型混为一谈。
    """
    if not text:
        return []
    stated: Dict[str, set] = defaultdict(set)
    for m in RE_ABSTRACT_ANOMER.finditer(text):
        code = canon_mono(m.group("sugar"))
        if code:
            stated[code].add("a" if m.group("anom") == "α" else "b")
    if not stated:
        return []
    seen: Dict[str, set] = defaultdict(set)
    for r in residues or []:
        code = canon_mono(_attr(r, "monosaccharide_name", "") or "")
        anom = _attr(r, "anomer")
        if code and anom in ("a", "b"):
            seen[code].add(anom)
    findings = []
    for code, st in stated.items():
        me = seen.get(code)
        if me and not (st & me):
            findings.append({
                "rule": "C5", "severity": "flagged",
                "message": (f"构型冲突：正文写 {sorted(st)}-{code}，"
                            f"但归属表中 {code} 残基实测均为 {sorted(me)}"
                            "（摘要与自表数据矛盾，需人工确认）"),
                "fields": {"metric": "anomer_conflict", "mono": code,
                           "stated": sorted(st), "observed": sorted(me)},
            })
    return findings


def check_methylation_counts(methylation: Dict[str, Any],
                             residues: Iterable[Any],
                             ratio: float = 2.5) -> List[Dict[str, Any]]:
    """C6：甲基化给出的连接类型数量 与 位移表残基类型数量 是否同一量级。

    两套数据本来不必逐条相等（甲基化覆盖全部残基、位移表只到可分辨的残基
    类型），但相差数倍通常意味着其中一套解析不全。
    """
    if not methylation:
        return []
    n_meth = len(methylation)
    n_res = len(list(residues or []))
    if not n_res:
        return []
    if max(n_meth, n_res) / max(1, min(n_meth, n_res)) >= ratio:
        return [{
            "rule": "C6", "severity": "warning",
            "message": (f"甲基化连接类型 {n_meth} 种 vs 位移归属残基 {n_res} 个，"
                        f"相差 ≥{ratio:g} 倍（可能有一套解析不全）"),
            "fields": {"metric": "methylation_vs_residue_counts",
                       "methylation_types": n_meth, "residue_entries": n_res},
        }]
    return []


def _get(obj: Any, key: str, default=None):
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def _residues_from(obj: Any) -> List[Any]:
    """从 core.GlycanRecord 或 dry-run dict 取残基列表。"""
    res = _get(obj, "residues") or []
    return list(res)


def _methylation_from(obj: Any) -> Dict[str, Any]:
    """从 poly_props.branching 提取甲基化连接类型（如 '4-GalA(p), t-GalA(p)'）。"""
    pp = _get(obj, "poly_props")
    if not pp:
        return {}
    branching = _get(pp, "branching") or ""
    types: Dict[str, Any] = {}
    # 注意：连接类型内部含逗号（如 "4,6-GalA(p)"），只能按分号/顿号切分；
    # 按逗号切会把 "4,6-GalA(p)" 拆成 "6-GalA(p)"（实测踩过）。
    # 分隔符是逗号，但连接类型内部也含逗号（"4,6-GalA(p)"），
    # 因此只在"逗号后紧跟字母"处切分（"4,6" 后面是数字，不会被切开）。
    for part in re.split(r",\s+(?=[A-Za-z0-9(])", branching):
        part = part.strip()
        if not part:
            continue
        m = re.match(r"^(?P<pos>\d+(?:\s*,\s*\d+)*|t)\s*-\s*(?P<sugar>[A-Za-z0-9]+)", part)
        if m:
            types[part] = m.group("sugar")
    return types


def cross_check(rec: Any,
                percentages: Optional[Dict[str, Any]] = None,
                statement_text: Optional[str] = None,
                tol_pct: float = 5.0) -> Dict[str, Any]:
    """对一条记录做全部交叉一致性校验。

    参数
    ----
    rec : ``core.GlycanRecord``（或等价的 dry-run dict）。从中读取
          residues / poly_props / molecular_formula / structure_level / composition。
    percentages : 单糖组成百分比；缺省从 ``poly_props.monosaccharide_ratio`` 取。
    statement_text : 摘要/正文陈述文本（用于 C5 构型冲突）。

    返回
    ----
    ``{"findings": [...], "flagged": bool, "n_flagged": int, "summary": str}``
    """
    residues = _residues_from(rec)
    pp = _get(rec, "poly_props")
    if percentages is None:
        percentages = (_get(pp, "monosaccharide_ratio") if pp else None) or {}
    methylation = _methylation_from(rec)

    findings: List[Dict[str, Any]] = []
    findings += check_composition_percent_sum(percentages, tol=tol_pct)
    findings += check_component_sets(percentages, residues, methylation=methylation)
    findings += check_formula_vs_residues(
        _get(rec, "molecular_formula"),
        _get(rec, "composition") or None,
        structure_level=_get(rec, "structure_level"))
    if statement_text:
        findings += check_anomer_conflict(statement_text, residues)
    findings += check_methylation_counts(methylation, residues)

    n_flagged = sum(1 for f in findings if f["severity"] == "flagged")
    return {
        "findings": findings,
        "flagged": n_flagged > 0,
        "n_flagged": n_flagged,
        "summary": "; ".join(f"[{f['rule']}] {f['message']}" for f in findings),
    }


def apply_to_record(rec: Any, findings: List[Dict[str, Any]],
                    prefix: str = "一致性") -> int:
    """把校验发现写入记录的 qc_notes，并在有 flagged 时置 qc_status='flagged'。

    返回写入的条数（便于测试断言）。
    """
    notes = getattr(rec, "qc_notes", None)
    if notes is None:
        return 0
    n = 0
    for f in findings:
        notes.append(f"{prefix}[{f['rule']}]: {f['message']}")
        n += 1
        if f["severity"] == "flagged":
            try:
                rec.qc_status = "flagged"
            except AttributeError:
                pass
    return n
