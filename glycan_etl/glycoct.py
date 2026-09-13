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
* ``composition_only``  —— 仅有残基组成，**不生成 GlycoCT**，改给组成式

下游可据此判断一条记录到底"能比对到什么程度"。

依赖
----
glypy 仅用于**校验**生成结果，缺失时不影响生成功能。
"""
from __future__ import annotations

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


def build_glycoct(residues: Sequence[Any],
                  sugar_type: str = "oligo",
                  chains: Optional[List[List[int]]] = None,
                  branch_links: Optional[List[Tuple[int, int, int]]] = None,
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

    返回
    ----
    dict(glycoct, level, composition, warnings, validated)
    """
    warnings: List[str] = []
    residues = list(residues)
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
