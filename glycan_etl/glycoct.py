# -*- coding: utf-8 -*-
"""标准 GlycoCT 编码生成与校验（P0-1）。

背景
----
早期版本在 core.py 里直接拼字符串产出所谓 "GlycoCT"，例如::

    RES 1b:a-lglcp-1:5|2:x

这不是合法 GlycoCT：缺少 RES/LIN 分段、basetype 误用三字母码
（glcp/glcpnac）、环信息写法错误。glypy、GlyTouCan 等工具都无法解析，
AI 平台也就无法用它做结构比对。

本模块按 GlycoCT 规范从「残基表 + 键连信息」生成标准编码：:

    RES
    1b:a-dglc-HEX-1:5
    2b:x-dglc-HEX-1:5
    LIN
    1:2o(4+1)1d

连接方向遵循 GlycoCT 惯例 ``<受体残基>o(<受体位点>+1)<供体残基>d``：
供体提供异头碳，受体提供羟基氧。糖醛酸写作 ``|6:a``、6-脱氧糖写作
``|6:d``；N-乙酰氨基糖按规范以独立取代基残基 ``s:n-acetyl`` 表达。

三级结构表达（关键设计）
------------------------
GlycoCT 只能表达**确定结构**。文献里大量多糖只给出甲基化/组成信息，
残基之间的连接顺序未知，此时硬生成 GlycoCT 等于伪造数据。因此本模块
区分三级：

* ``complete``          —— 单糖 / 寡糖，连接明确，生成完整 GlycoCT
* ``repeat_unit``       —— 单一重复单元多糖，生成一个重复单元的 GlycoCT
* ``domain_only``       —— 只给出**域级骨架**（如"主链为 HG 域 + 少量带侧链的 RG-I 域"），
                          残基到残基的顺序未知，**不生成 GlycoCT**，改给 domain_architecture
* ``composition_only``  —— 仅有残基组成，**不生成 GlycoCT**，改给组成式

`domain_only` 是 2025-09 新增的一档：很多多糖文献（尤其果胶类）只在正文写一句
"主链由大量 HG 域与少量带侧链的 RG-I 域组成"，既没有逐残基连接式，也不只是组成
百分比。硬塞进 composition_only 会丢掉"主链是什么、侧链挂在哪个域上"这一层信息；
硬生成 GlycoCT 又是伪造。因此单独设一档，携带结构化域架构。

下游可据此判断一条记录到底"能比对到什么程度"。

依赖
----
glypy 仅用于**校验**生成结果，缺失时不影响生成功能。
"""
from __future__ import annotations

import re
from collections import Counter
from typing import Any, Dict, List, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# 单糖定义表：code -> (绝对构型, basetype, superclass, 环, 修饰段)
# ---------------------------------------------------------------------------
# 构型必须逐个显式给定：Fuc / Rha / IdoA 是 L 构型，不能用母体糖（D 构型）
# 推导，否则会上游产生错误的立体化学。
# 局限：文献只写 "Ara" / "Gal" 而没标 D/L 时，这里取最常见的天然构型，
#       属于已知近似（可后续从原文标题的 "L-Araf" 等写法补全）。
MONO_DEF: Dict[str, Tuple[str, str, str, str, str]] = {
    # 己醛糖
    "Glc": ("d", "glc", "HEX", "1:5", ""),
    "Gal": ("d", "gal", "HEX", "1:5", ""),
    "Man": ("d", "man", "HEX", "1:5", ""),
    "All": ("d", "all", "HEX", "1:5", ""),
    "Alt": ("d", "alt", "HEX", "1:5", ""),
    "Gul": ("d", "gul", "HEX", "1:5", ""),
    "Ido": ("l", "ido", "HEX", "1:5", ""),
    "Tal": ("d", "tal", "HEX", "1:5", ""),
    # 戊醛糖
    "Xyl": ("d", "xyl", "PEN", "1:4", ""),
    "Ara": ("l", "ara", "PEN", "1:4", ""),
    "Rib": ("d", "rib", "PEN", "1:4", ""),
    "Lyx": ("d", "lyx", "PEN", "1:4", ""),
    # 酮糖（glypy 的 stem 字典不含酮糖，生成后无法用 glypy 校验）
    "Fru": ("d", "fru", "HEX", "2:5", ""),
    # 糖醛酸：C6 为羧酸
    "GlcA": ("d", "glc", "HEX", "1:5", "|6:a"),
    "GalA": ("d", "gal", "HEX", "1:5", "|6:a"),
    "ManA": ("d", "man", "HEX", "1:5", "|6:a"),
    "IdoA": ("l", "ido", "HEX", "1:5", "|6:a"),
    "GulA": ("d", "gul", "HEX", "1:5", "|6:a"),
    # 6-脱氧糖：C6 为甲基（Fuc = 6-脱氧-L-半乳糖，Rha = 6-脱氧-L-甘露糖）
    "Fuc": ("l", "gal", "HEX", "1:5", "|6:d"),
    "Rha": ("l", "man", "HEX", "1:5", "|6:d"),
    "Qui": ("d", "glc", "HEX", "1:5", "|6:d"),
}

# N-乙酰氨基糖：以独立取代基残基表达（GlycoCT 规范写法，glypy 可解析）
NACETYL_DEF: Dict[str, Tuple[str, str, str]] = {
    "GlcNAc": ("d", "glc", "1:5"),
    "GalNAc": ("d", "gal", "1:5"),
    "ManNAc": ("d", "man", "1:5"),
}

# 异头构型：对外统一用 a/b/x
ANOMER_MAP: Dict[Optional[str], str] = {
    "a": "a", "b": "b", "x": "x",
    "alpha": "a", "beta": "b", "unknown": "x", None: "x", "": "x",
}

# 常见单糖的"异头碳位置"（醛糖 C1，酮糖 C2）
ANOMERIC_CARBON = {"Fru": 2, "Sorb": 2}


def _anomer_char(anomer: Optional[str]) -> str:
    return ANOMER_MAP.get((anomer or "").strip().lower(), "x")


def _residue_line(idx: int, name: str, anomer: Optional[str]) -> Optional[str]:
    """生成一行 RES 内容；未收录的单糖返回 None。"""
    if name in NACETYL_DEF:
        config, stem, ring = NACETYL_DEF[name]
        return f"{idx}b:{_anomer_char(anomer)}-{config}{stem}-HEX-{ring}"
    if name in MONO_DEF:
        config, stem, sup, ring, mod = MONO_DEF[name]
        return f"{idx}b:{_anomer_char(anomer)}-{config}{stem}-{sup}-{ring}{mod}"
    return None


def glycan_composition(residues: Sequence[Any]) -> str:
    """残基组成式，如 ``GalA6,Gal3,Ara3,Rha2,GlcA``（按数量降序）。"""
    counts: Counter = Counter()
    for r in residues:
        name = (getattr(r, "monosaccharide_name", "") or "").strip()
        if name:
            counts[name] += 1
    if not counts:
        return ""
    return ",".join(f"{k}{v}" if v > 1 else k
                    for k, v in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))


def _composition_only(residues: Sequence[Any], warnings: List[str]) -> Dict[str, Any]:
    return {
        "glycoct": None,
        "level": "composition_only",
        "composition": glycan_composition(residues),
        "warnings": warnings,
        "validated": False,
    }


# ---------------------------------------------------------------------------
# 域级骨架（domain_only）：从正文结论句抽取"主链是哪个域、侧链挂在哪"
# ---------------------------------------------------------------------------
# 果胶类多糖的标准域名词（HG / RG-I / RG-II / XGA 等），文献里几乎只用这些写法。
DOMAIN_ALIASES: Dict[str, str] = {
    "hg": "HG", "homogalacturonan": "HG", "homogalacturonans": "HG",
    "rg-i": "RG-I", "rgi": "RG-I", "rg-1": "RG-I", "rg1": "RG-I",
    "rhamnogalacturonan-i": "RG-I", "rhamnogalacturonan i": "RG-I",
    "rg-ii": "RG-II", "rgii": "RG-II", "rg-2": "RG-II", "rg2": "RG-II",
    "rhamnogalacturonan-ii": "RG-II", "rhamnogalacturonan ii": "RG-II",
    "xga": "XGA", "xylogalacturonan": "XGA",
    "aga": "AGA", "arabinogalactan": "AG", "ag": "AG",
    "arabinogalactan-i": "AG-I", "arabinogalactan-ii": "AG-II",
    "arabinan": "arabinan", "arabinan domain": "arabinan",
    "galactan": "galactan", "galacturonan": "galacturonan",
}

# "composed of ... HG domains ..." 类结论句：必须有域名词才算命中
RE_DOMAIN_SENTENCE = re.compile(
    r"([^.]{0,160}?\b(?:composed\s+of|consists?\s+of|comprising|"
    r"mainly\s+composed\s+of|dominated\s+by)[^.]{0,200}?"
    r"\b(?:hg|homogalacturonan|rg[\s\-]?(?:i{1,2}|[12])|xga|"
    r"rhamnogalacturonan[\s\-]?(?:i{1,2}|[12])|arabinogalactan|"
    r"arabinan|galactan)\b[^.]{0,200})", re.I)

RE_DOMAIN_TOKEN = re.compile(
    r"\b(?P<qty>large|small|major|minor|a\s+few|several|significant)?\s*"
    r"(?P<name>homogalacturonan|rhamnogalacturonan[\s\-]?(?:i{1,2}|[12])|"
    r"rg[\s\-]?(?:i{1,2}|[12])|hg|xga|arabinogalactan|arabinan|galactan)\b"
    r"(?P<domains>\s+domains?)?", re.I)

RE_QUANTITY = re.compile(
    r"\b(large\s+(?:number|amount)\s+of|small\s+(?:number|amount)\s+of|"
    r"major|minor|a\s+few|several|significant\s+amount\s+of)\b", re.I)


def _canon_domain(name: str) -> str:
    key = re.sub(r"\s+", " ", (name or "").strip().lower())
    key = key.replace("–", "-").replace("—", "-")
    if key in DOMAIN_ALIASES:
        return DOMAIN_ALIASES[key]
    # "rhamnogalacturonan i" / "rg i" 这类带空格或罗马数字的写法
    key2 = re.sub(r"[\s\-]+", "-", key)
    if key2 in DOMAIN_ALIASES:
        return DOMAIN_ALIASES[key2]
    key3 = key2.replace("-", "")
    if key3 in DOMAIN_ALIASES:
        return DOMAIN_ALIASES[key3]
    return name.strip()


def extract_domain_architecture(text: str) -> Dict[str, Any]:
    """从正文抽取**域级骨架**结论，如：

        "HP is mainly composed of a large number of HG domains and a
         small number of RG-I domains with side chains."

    → ``{"domains": [{"name": "HG", "quantity": "large", "side_chains": False},
                     {"name": "RG-I", "quantity": "small", "side_chains": True}],
         "evidence": "<原句>", "source": "text_conclusion"}``

    为什么需要这一档：很多多糖文献（尤其果胶）既不给逐残基连接式，也不只是
    组成百分比，而是给一句域级结论。这类句子不含 '→'，因此
    :func:`table_parser.extract_linkage_sequence` 的触发词抓不到，导致整条
    结构信息被丢弃、记录降级为 composition_only。
    """
    if not text:
        return {"domains": []}
    t = re.sub(r"\s+", " ", text)
    m = RE_DOMAIN_SENTENCE.search(t)
    if not m:
        return {"domains": []}
    sentence = m.group(1).strip()
    domains: List[Dict[str, Any]] = []
    seen: set = set()
    for dm in RE_DOMAIN_TOKEN.finditer(sentence):
        raw = dm.group("name")
        # 裸 "rg" 不可信（可能是别的缩写），要求带罗马数字/阿拉伯数字
        if re.fullmatch(r"rg", raw.strip(), re.I):
            continue
        name = _canon_domain(raw)
        if not name:
            continue
        # 左右各取一段用于判定数量词与"with side chains"
        lo, hi = dm.start(), dm.end()
        left = sentence[max(0, lo - 40):lo]
        right = sentence[hi:hi + 40]
        qty = None
        qm = RE_QUANTITY.search(left)
        if qm:
            q = qm.group(1).lower()
            qty = ("large" if q.startswith(("large", "major", "significant"))
                   else "small" if q.startswith(("small", "minor"))
                   else "few")
        if qty is None:
            qm2 = RE_QUANTITY.search(right)
            if qm2:
                q = qm2.group(1).lower()
                qty = ("large" if q.startswith(("large", "major", "significant"))
                       else "small" if q.startswith(("small", "minor")) else "few")
            elif dm.group("qty"):
                q = dm.group("qty").lower()
                qty = ("large" if q.startswith(("large", "major", "significant"))
                       else "small" if q.startswith(("small", "minor")) else "few")
        side = bool(re.search(r"with\s+(?:side|branch)", right, re.I))
        key = (name, qty, side)
        if key in seen:
            continue
        seen.add(key)
        domains.append({"name": name, "quantity": qty, "side_chains": side})
    if not domains:
        return {"domains": []}
    return {"domains": domains, "evidence": sentence, "source": "text_conclusion"}



def build_glycoct(residues: Sequence[Any],
                  sugar_type: str = "oligo",
                  chains: Optional[List[List[int]]] = None,
                  branch_links: Optional[List[Tuple[int, int, int]]] = None,
                  domain_architecture: Optional[Dict[str, Any]] = None,
                  ) -> Dict[str, Any]:
    """从残基表生成标准 GlycoCT。

    参数
    ----
    residues : 具备 ``residue_seq`` / ``monosaccharide_name`` / ``anomer`` /
        ``parent_carbon`` 属性的对象序列（即 core.Residue）。
    sugar_type : ``mono`` / ``oligo`` / ``poly``，决定结构表达等级。
    chains : 可选的链定义，每个元素是**残基序号**列表（按 非还原端→还原端
        顺序）。缺省时把所有残基视为一条链。多条链用于主链 + 支链。
    branch_links : 支链挂接，元素为 ``(受体残基序号, 受体位点, 供体残基序号)``，
        表示供体的异头碳连到受体的 ``O<位点>``。
    domain_architecture : 见 :func:`extract_domain_architecture`。给出**域级骨架**
        （如 HG 主链 + 带侧链的 RG-I）时，残基间顺序通常仍未知，此时返回
        ``level='domain_only'`` 且 ``glycoct=None`` —— 既不伪造结构，也不把
        域级信息压缩成纯组成式。

    返回
    ----
    dict(glycoct, level, composition, warnings, validated)
    """
    warnings: List[str] = []
    residues = list(residues)
    if domain_architecture and domain_architecture.get("domains"):
        names = "、".join(d["name"] for d in domain_architecture["domains"])
        warnings.append(
            f"仅到域级骨架（{names}），残基间连接顺序未知，不生成 GlycoCT")
        return {
            "glycoct": None,
            "level": "domain_only",
            "composition": glycan_composition(residues),
            "domain_architecture": domain_architecture,
            "warnings": warnings,
            "validated": False,
        }
    if not residues:
        return _composition_only(residues, ["无残基信息，无法生成结构编码"])

    # 单残基重复单元：残基本身既是供体又是受体。若能从
    # parent_carbon(供体异头碳) + linkage_branch(受体羟基位点) 解出连接，
    # 就展开成"一个拷贝"的重复单元，把连接信息保留下来；否则只记组成。
    if sugar_type == "poly" and len(residues) == 1:
        import copy
        r0 = residues[0]
        acceptor_pos = getattr(r0, "linkage_branch", None)
        if acceptor_pos:
            twin = copy.copy(r0)
            twin.residue_seq = 2
            twin.parent_carbon = acceptor_pos     # 副本作为受体
            twin.is_reducing_end = False
            twin.linkage_branch = 0
            residues = [r0, twin]
        else:
            warnings.append(
                "单残基重复单元缺少受体位点(linkage_branch)，仅记录组成")
            return _composition_only(residues, warnings)

    res_lines: List[str] = []
    sub_lines: List[str] = []
    link_lines: List[str] = []
    gid: Dict[int, int] = {}          # residue_seq -> GlycoCT 残基编号
    next_id = 1

    # ---- RES 段 ----
    for r in residues:
        name = (getattr(r, "monosaccharide_name", "") or "").strip()
        line = _residue_line(next_id, name, getattr(r, "anomer", None))
        if line is None:
            warnings.append(f"未收录的单糖 {name!r}，无法生成 GlycoCT")
            return _composition_only(residues, warnings)
        gid[getattr(r, "residue_seq", next_id)] = next_id
        res_lines.append(line)
        next_id += 1
        # N-乙酰氨基糖需追加独立取代基残基 + 其连接
        if name in NACETYL_DEF:
            sub_lines.append(f"{next_id}s:n-acetyl")
            link_lines.append(f"{len(link_lines) + 1}:{gid[getattr(r, 'residue_seq', next_id)]}d(2+1){next_id}n")
            next_id += 1

    # ---- LIN 段：供体(异头碳) -> 受体(parent_carbon) ----
    by_seq = {getattr(r, "residue_seq", i + 1): r for i, r in enumerate(residues)}
    # 链定义：缺省为单链（所有残基按给定顺序）；多条链用于主链 + 支链
    chain_list: List[List[int]] = chains if chains else [
        [getattr(r, "residue_seq", i + 1) for i, r in enumerate(residues)]]

    chain_ok = True
    for chain in chain_list:
        for a_seq, b_seq in zip(chain, chain[1:]):
            acceptor = by_seq.get(b_seq)
            pos = getattr(acceptor, "parent_carbon", None) if acceptor else None
            if acceptor is None or pos is None:
                chain_ok = False
                warnings.append(f"残基 {b_seq} 缺少连接位点，无法确定连接顺序")
                break
            link_lines.append(
                f"{len(link_lines) + 1}:{gid[b_seq]}o({pos}+1){gid[a_seq]}d")
        if not chain_ok:
            break

    # 支链挂接：供体的异头碳 -> 受体的 O<位点>
    if chain_ok and branch_links:
        for acc_seq, acc_pos, don_seq in branch_links:
            if acc_seq not in gid or don_seq not in gid or not acc_pos:
                chain_ok = False
                warnings.append(
                    f"支链挂接 ({acc_seq}, {acc_pos}, {don_seq}) 无法解析")
                break
            link_lines.append(
                f"{len(link_lines) + 1}:{gid[acc_seq]}o({acc_pos}+1){gid[don_seq]}d")

    if not chain_ok:
        return _composition_only(residues, warnings)

    # ---- 组装 ----
    parts = ["RES"] + res_lines + sub_lines
    if link_lines:
        parts += ["LIN"] + link_lines
    text = "\n".join(parts)

    level = "complete" if sugar_type in ("mono", "oligo") else "repeat_unit"
    result = {
        "glycoct": text,
        "level": level,
        "composition": glycan_composition(residues),
        "warnings": warnings,
        "validated": False,
    }
    ok, msg = validate_glycoct(text)
    result["validated"] = bool(ok)
    if ok is False:
        result["warnings"].append(f"GlycoCT 校验失败: {msg}")
    elif ok is None:
        result["warnings"].append(msg)
    return result


def validate_glycoct(text: str) -> Tuple[Optional[bool], str]:
    """用 glypy 校验 GlycoCT 合法性。

    返回 ``(True, 'ok')`` / ``(False, 原因)`` / ``(None, 原因)``（无 glypy）。
    """
    if not text:
        return False, "空编码"
    try:
        from glypy.io import glycoct as _glycoct  # type: ignore
    except ImportError:
        return None, "未安装 glypy，跳过结构校验"
    try:
        _glycoct.loads(text)
        return True, "ok"
    except Exception as exc:  # noqa: BLE001
        return False, f"{type(exc).__name__}: {exc}"
