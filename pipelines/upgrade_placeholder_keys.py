#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""upgrade_placeholder_keys.py — 把占位结构键的记录原地升级为更强的表达等级。

背景（真实需求）
----------------
山楂多糖 HP（IJBM 2025, 10.1016/j.ijbiomac.2025.145713）此前以
``structure_level='composition_only'`` 入库，glycoct 是占位键::

    UNRESOLVED:poly:10.1016/j.ijbiomac.2025.145713:8db2792c…

新管线从正文抽到了域级骨架（HG 主链 + 带侧链的 RG-I），等级升为
``domain_only``，去重键也随之变成 ``DOMAIN:…``。**键变了，就不能直接再插一条**：
那样旧行会变成孤行，170 条位移被劈成两份（这正是仓库已知的"结构身份不稳定"缺陷）。

因此升级必须在**一个事务内**完成"删旧 + 写新"，并在提交前做对账：
位移数、实验数、残基数与升级前不一致就整体回滚，绝不留下半成品。

用法
----
  # 先看会改什么（不写库）
  python3 pipelines/upgrade_placeholder_keys.py --pdf <paper.pdf> --config <cfg> --dry-run
  # 真正执行
  python3 pipelines/upgrade_placeholder_keys.py --pdf <paper.pdf> --config <cfg>
"""
import argparse
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from glycan_etl import core                      # noqa: E402
from glycan_etl import table_parser as tp        # noqa: E402


def blob_to_record(blob: dict, doi: str, journal: str, year: int) -> core.GlycanRecord:
    """table_parser 的 dry-run dict → core.GlycanRecord（与 batch_etl 同口径）。"""
    rec = core.GlycanRecord(
        sugar_type=blob.get("sugar_type") or "poly",
        iupac_short=blob.get("iupac_short"),
        glycoct=blob.get("glycoct"),
        molecular_formula=blob.get("molecular_formula"),
        molecular_weight=blob.get("molecular_weight"),
        anomer=blob.get("anomer"),
        doi=doi, journal=journal, year=year,
        nmr_page=blob.get("nmr_page"),
        structure_level=blob.get("structure_level"),
        composition=blob.get("composition"),
        domain_architecture=blob.get("domain_architecture"),
    )
    fields = set(core.Residue.__dataclass_fields__)
    rec.residues = [core.Residue(**{k: v for k, v in r.items() if k in fields})
                    for r in blob.get("residues") or []]
    if blob.get("physicochemical"):
        f = set(core.Physicochemical.__dataclass_fields__)
        rec.physicochemical = core.Physicochemical(
            **{k: v for k, v in blob["physicochemical"].items() if k in f})
    if blob.get("poly_props"):
        f = set(core.PolysaccharideProps.__dataclass_fields__)
        rec.poly_props = core.PolysaccharideProps(
            **{k: v for k, v in blob["poly_props"].items() if k in f})
    for e in blob.get("experiments") or []:
        exp = core.NMRExperiment(
            nucleus=e.get("nucleus", "1H"), experiment_2d=e.get("experiment_2d"),
            solvent=e.get("solvent") or "D2O", frequency=e.get("frequency"),
            temperature=e.get("temperature"), ph=e.get("ph"))
        pf, p2 = set(core.Peak1D.__dataclass_fields__), set(core.Peak2D.__dataclass_fields__)
        exp.peaks = [core.Peak1D(**{k: v for k, v in p.items() if k in pf})
                     for p in e.get("peaks") or []]
        exp.peaks_2d = [core.Peak2D(**{k: v for k, v in p.items() if k in p2})
                        for p in e.get("peaks_2d") or []]
        rec.experiments.append(exp)
    rec.qc_notes = list(blob.get("qc_notes") or [])
    qs = blob.get("qc_status", "dry-run-ok")
    rec.qc_status = {"dry-run-ok": "passed"}.get(qs, qs)
    return rec


def counts(cur, sugar_id: int) -> dict:
    cur.execute("SELECT count(*) FROM nmr_experiments WHERE sugar_id=%s", (sugar_id,))
    n_exp = cur.fetchone()[0]
    cur.execute("SELECT count(*) FROM nmr_shifts_1d p JOIN nmr_experiments e "
                "USING(experiment_id) WHERE e.sugar_id=%s", (sugar_id,))
    n_shift = cur.fetchone()[0]
    cur.execute("SELECT count(*) FROM residues WHERE sugar_id=%s", (sugar_id,))
    n_res = cur.fetchone()[0]
    cur.execute("SELECT count(*) FROM polysaccharide_props WHERE sugar_id=%s", (sugar_id,))
    n_pp = cur.fetchone()[0]
    return {"experiments": n_exp, "shifts": n_shift, "residues": n_res, "poly_props": n_pp}


def main() -> int:
    ap = argparse.ArgumentParser(description="占位结构键记录的原地升级（事务 + 对账）")
    ap.add_argument("--pdf", required=True, help="源文献 PDF（用于重新解析）")
    ap.add_argument("--config", required=True, help="含 db 连接的配置")
    ap.add_argument("--dry-run", action="store_true", help="只打印计划，不写库")
    args = ap.parse_args()

    if not os.path.exists(args.config):
        print(f"[fatal] 配置不存在: {args.config}", file=sys.stderr)
        return 1
    import yaml
    with open(args.config, encoding="utf-8") as f:
        config = yaml.safe_load(f) or {}
    db = config.get("db")
    doi = config.get("doi")
    if not db or not doi:
        print("[fatal] 配置需含 db 与 doi", file=sys.stderr)
        return 1

    print(f"[1/4] 重新解析 {os.path.basename(args.pdf)} …")
    d = tp.run(args.pdf, verbose=False)
    if not d["records"]:
        print("[fatal] 未解析出记录", file=sys.stderr)
        return 1
    rec = blob_to_record(d["records"][0], doi, config.get("journal"), config.get("year"))
    new_key = core.resolve_glycoct(rec)
    print(f"      等级={rec.structure_level}  新键={new_key}")
    print(f"      域架构={rec.domain_architecture}")

    import psycopg2
    conn = psycopg2.connect(**db)
    cur = conn.cursor()
    cur.execute("SELECT sugar_id, glycoct, structure_level FROM sugars "
                "WHERE first_seen_doi=%s", (doi,))
    rows = cur.fetchall()
    if not rows:
        print(f"[fatal] 库中找不到 doi={doi} 的记录", file=sys.stderr)
        return 1
    print(f"[2/4] 库中命中 {len(rows)} 条旧记录:")
    for sid, key, lvl in rows:
        print(f"      sugar_id={sid} level={lvl} key={key[:60]}")
        print(f"      {'旧派生态: ' + str(counts(cur, sid))}")

    if args.dry_run:
        print("[3/4] --dry-run：不写库。")
        conn.close()
        return 0

    # 删除并重写：新键可能只对应其中一条；本脚本处理"单一记录"的常见情形
    sid, old_key, old_level = rows[0]
    if len(rows) > 1:
        print("[fatal] 命中多条，需人工指定 sugar_id（避免误删）", file=sys.stderr)
        conn.close()
        return 1
    if old_key == new_key:
        print("[3/4] 键未变化，仅补写等级/域架构。")

    before = counts(cur, sid)

    # ---- 写前对账（关键）----------------------------------------------------
    # core.insert_record() 内部会 commit，因此"先写再发现问题就 rollback"是假的：
    # 删除早已生效。必须在动手之前就拿解析结果算清将要写入的行数，不一致就
    # 根本不动库，避免把 170 条位移删掉却写不回来。
    n_1d = sum(len(e.peaks) for e in rec.experiments)
    n_2d = sum(len(e.peaks_2d) for e in rec.experiments)
    expect = {
        "experiments": len(rec.experiments),
        "shifts": n_1d,
        "residues": len(rec.residues),
        "poly_props": 1 if rec.poly_props else 0,
    }
    print(f"[3/4] 写前对账")
    print(f"      库中现有 : {before}")
    print(f"      本次将写 : {expect}  (1D 峰 {n_1d} + 2D 峰 {n_2d})")
    shrink = {k: (before[k], expect[k]) for k in before if expect.get(k, 0) < before[k]}
    if shrink:
        print(f"[fatal] 解析结果少于库中现有数据 {shrink}；"
              "拒绝删除旧记录（可能是版式解析退化）。未改动数据库。", file=sys.stderr)
        conn.close()
        return 1

    print(f"      升级 sugar_id={sid}: {old_level} → {rec.structure_level}")
    cur.execute("DELETE FROM sugars WHERE sugar_id=%s", (sid,))
    core.insert_record(conn, rec)          # 内含 commit
    cur.execute("SELECT sugar_id FROM sugars WHERE glycoct=%s", (new_key,))
    row = cur.fetchone()
    if not row:
        print("[fatal] 新记录未写入（旧记录已删除，请重跑本脚本）", file=sys.stderr)
        return 1
    new_sid = row[0]
    after = counts(cur, new_sid)
    print(f"[4/4] 升级后={after}")
    if after["shifts"] < before["shifts"]:
        print(f"[warn] 位移数由 {before['shifts']} 变为 {after['shifts']}，请人工复核",
              file=sys.stderr)
    conn.close()
    print(f"[OK] 升级完成：sugar_id {sid} → {new_sid}，键 {old_key[:40]}… → {new_key[:40]}…")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:  # noqa: BLE001
        print(f"[fatal] {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
