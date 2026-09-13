#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""backfill_glycoct.py — 回填/修正 sugars 表的结构编码（P0-1）。

解决的问题
----------
早期版本写入的是手拼伪编码（如 ``RES 1b:a-lglcp-1:5|2:x``），不是合法
GlycoCT，外部工具无法解析。本脚本按现有 residues 重新生成标准编码，
并补齐 structure_level / composition。

处理规则
--------
* 能生成确定结构  -> 用标准 GlycoCT 覆盖旧编码（同一结构的多条记录若因此
  撞上 UNIQUE 约束，会被报告出来，默认**不自动合并**，需人工确认）。
* 只有组成信息    -> 保留原占位编码（sugars.glycoct 为 NOT NULL UNIQUE），
  但写入 structure_level='composition_only' 与 composition，并清空误导性的
  伪结构编码（替换为可识别的 COMPOSITION 占位串）。

用法
----
  # 先看会改什么（推荐）
  python3 pipelines/backfill_glycoct.py --config pipelines/config.local.yaml --dry-run
  # 真正写入
  python3 pipelines/backfill_glycoct.py --config pipelines/config.local.yaml
"""
import argparse
import hashlib
import os
import re
import sys
from collections import defaultdict

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from glycan_etl import glycoct as gx_lib          # noqa: E402
from glycan_etl import core as core_lib           # noqa: E402
from glycan_etl.core import GlycanRecord, Residue, resolve_glycoct   # noqa: E402


def guess_mono(name: str):
    """从 iupac_short 猜单糖码。

    先走标题匹配；再处理 "GlcpNAc" 这类把环形式标记 p/f 插在糖名中间、
    导致 ``\\bGlcNAc\\b`` 匹配不上的写法；最后退化为包含匹配。
    """
    code = core_lib.parse_mono_title(name)
    if code:
        return code
    norm = re.sub(r"(?<=[A-Za-z])[pf](?=[A-Z])", "", name)
    code = core_lib.parse_mono_title(norm)
    if code:
        return code
    low = norm.lower()
    for c in sorted(core_lib.MONOSACCHARIDES, key=len, reverse=True):
        if c.lower() in low:
            return c
    return None


def residues_for(sugar_id: int, iupac_short: str, sugar_type: str, res_map: dict) -> list:
    """取残基列表；老记录可能根本没有 residues 行，则从 iupac_short 重建。

    早期 v001 种子数据只写了 sugars / nmr_experiments / nmr_shifts_1d，
    residues 表是 v002 才引入的，因此这些记录没有残基信息。直接回填会把
    它们降级成 "COMPOSITION:unknown"，比原编码信息更少；这里改为用解析器
    从 iupac_short 还原残基。
    """
    residues = res_map.get(sugar_id, [])
    if residues:
        return residues
    if not iupac_short:
        return []
    # 单糖：直接从名称还原（build_residues 对 "α-D-Glcp" 这类带构型前缀的
    # 写法会因 \\b 词边界匹配不上而返回空列表）
    if (sugar_type or "mono") == "mono":
        code = guess_mono(iupac_short)
        if code:
            anom = core_lib.parse_mono_anomer(iupac_short) or "unknown"
            return [Residue(1, code, ring_form="p", anomer=anom, is_reducing_end=True)]
    try:
        residues = core_lib.build_residues(iupac_short, sugar_type or "mono")
    except Exception:  # noqa: BLE001
        return []
    # 标题里的 α/β（如 "α-D-Glcp"）比糖名本身更可靠，回填到残基上
    anom = core_lib.parse_mono_anomer(iupac_short)
    if anom:
        for r in residues:
            if (r.anomer or "").strip().lower() in ("", "unknown", "x"):
                r.anomer = anom
    return residues


def looks_like_legacy_glycoct(text: str) -> bool:
    """旧版伪编码的判别（缺 RES/LIN 分段，或 basetype 用了三字母码）。"""
    if not text:
        return False
    if text.startswith("UNRESOLVED:") or text.startswith("COMPOSITION:"):
        return False
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if "RES" not in lines:          # 标准 GlycoCT 必有 RES 分段
        return True
    return False


def build_plan(cur) -> list:
    """读出全部 sugar + residues, 生成逐条更新计划。"""
    cur.execute(
        "SELECT sugar_id, sugar_type, iupac_short, glycoct, first_seen_doi "
        "FROM sugars ORDER BY sugar_id")
    sugars = cur.fetchall()
    cur.execute(
        "SELECT sugar_id, residue_seq, monosaccharide_name, ring_form, anomer, "
        "is_reducing_end, parent_carbon, linkage_branch FROM residues "
        "ORDER BY sugar_id, residue_seq")
    res_map = defaultdict(list)
    for row in cur.fetchall():
        res_map[row[0]].append(Residue(
            residue_seq=row[1], monosaccharide_name=row[2] or "", ring_form=row[3] or "p",
            anomer=row[4] or "unknown", is_reducing_end=bool(row[5]),
            parent_carbon=row[6], linkage_branch=row[7] or 0))

    plan = []
    for sid, stype, name, old, doi in sugars:
        residues = residues_for(sid, name, stype, res_map)
        built = gx_lib.build_glycoct(residues, stype or "poly")
        new_glycoct = built["glycoct"]
        level, comp = built["level"], built["composition"]
        if new_glycoct is None:
            # 组成型结构：沿用 ETL 的占位键规则（结构指纹 + DOI），
            # 保证回填后再跑 ETL 仍命中同一条记录，不会因键变化而重复建库。
            new_glycoct = resolve_glycoct(GlycanRecord(
                sugar_type=stype or "poly", iupac_short=name, doi=doi,
                residues=residues, glycoct=None))
        changed = (new_glycoct != old) or (level is not None)
        plan.append({
            "sugar_id": sid, "name": name, "sugar_type": stype,
            "old": old, "new": new_glycoct, "level": level, "composition": comp,
            "n_residues": len(residues),
            "legacy": looks_like_legacy_glycoct(old),
            "changed": changed,
        })
    return plan


def main() -> int:
    ap = argparse.ArgumentParser(description="回填标准 GlycoCT 与结构表达等级")
    ap.add_argument("--config", default="pipelines/config.local.yaml",
                    help="含 db 连接的配置文件")
    ap.add_argument("--dry-run", action="store_true", help="只打印计划，不写库")
    args = ap.parse_args()

    if not os.path.exists(args.config):
        print(f"[fatal] 配置文件不存在: {args.config}", file=sys.stderr)
        return 1
    import yaml
    with open(args.config, encoding="utf-8") as f:
        config = yaml.safe_load(f) or {}
    db = config.get("db")
    if not db:
        print("[fatal] 配置缺少 db 连接信息", file=sys.stderr)
        return 1

    try:
        import psycopg2
    except ImportError:
        print("[fatal] 需要 psycopg2-binary", file=sys.stderr)
        return 1

    conn = psycopg2.connect(**db)
    cur = conn.cursor()
    plan = build_plan(cur)

    n_legacy = sum(1 for p in plan if p["legacy"])
    n_comp = sum(1 for p in plan if p["level"] == "composition_only")
    print(f"共 {len(plan)} 条结构记录；其中旧版伪编码 {n_legacy} 条，"
          f"仅组成型 {n_comp} 条\n")

    for p in plan:
        flag = "LEGACY" if p["legacy"] else "      "
        mark = "改" if p["changed"] else "留"
        print(f"  [{mark}] {flag} #{p['sugar_id']:<3} {str(p['name'])[:34]:34s} "
              f"level={p['level']:<16s} res={p['n_residues']}")
        if p["new"] != p["old"]:
            print(f"        旧: {(p['old'] or '')[:70]}")
            print(f"        新: {p['new'][:70].replace(chr(10), ' | ')}")

    if args.dry_run:
        print("\n[dry-run] 未写库。去掉 --dry-run 即执行。")
        conn.close()
        return 0

    # 唯一键冲突预检：新编码若与别的记录重复，说明是同一结构，需人工合并
    seen, conflicts = {}, []
    for p in plan:
        key = p["new"]
        if key in seen and seen[key] != p["sugar_id"]:
            conflicts.append((seen[key], p["sugar_id"], key))
        else:
            seen[key] = p["sugar_id"]
    if conflicts:
        print(f"\n[!] 发现 {len(conflicts)} 组编码冲突（同一结构多条记录），"
              "已跳过这些记录，请人工确认合并策略：")
        for a, b, key in conflicts[:10]:
            print(f"    sugar_id {a} / {b} -> {key[:60]}")
    conflict_ids = {sid for c in conflicts for sid in c[:2]}

    n_done = 0
    for p in plan:
        if p["sugar_id"] in conflict_ids or not p["changed"]:
            continue
        cur.execute(
            "UPDATE sugars SET glycoct=%s, glycoct_hash=encode(digest(%s,'sha256'),'hex'), "
            "structure_level=COALESCE(%s, structure_level), "
            "composition=COALESCE(%s, composition), updated_at=now() "
            "WHERE sugar_id=%s",
            (p["new"], p["new"], p["level"], p["composition"], p["sugar_id"]),
        )
        n_done += 1
    conn.commit()
    cur.close()
    conn.close()
    print(f"\n完成: 更新 {n_done} 条，跳过 {len(conflict_ids)} 条冲突记录。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
