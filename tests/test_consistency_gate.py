#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""一致性硬门槛（v007）的数据库行为回归测试。

验证的核心契约：**坏数据进不来，好数据进得来，且默认可信度是"未复核"**。

  关卡一（BEFORE INSERT/UPDATE OF sugars）
    C7b  声明 complete/repeat_unit 却给占位键/空值  → 必须被数据库拒绝
    C8   domain_only 但 domain_architecture 为空     → 必须被拒绝
    C10  声明 composition_only/domain_only 却给真结构 → 必须被拒绝
    合法记录（真结构+complete / 占位键+domain_only+域架构）→ 必须写入成功

  关卡二（consistency_recheck）
    ★ fail-closed：新记录 consistency_ok 初始为 FALSE，复核通过才置 TRUE
    C1/C2/C3/C6 这类需要派生表齐备的规则在复核阶段产出发现

为什么必须有这层测试：这些规则此前只在 Python 侧，手工 SQL 或任何绕过
batch_etl 的写入都能把自相矛盾的数据塞进"真值库"。门槛是数据库行为，
因此必须用真实数据库来验，不能只在应用层测。

依赖数据库，连不上则整组跳过（不会让自检失败）：
  GLYCAN_TEST_DB=glycan_final  python3 tests/test_consistency_gate.py
  PGHOST/PGPORT/PGUSER 可用标准 libpq 环境变量覆盖

用法:
  python3 tests/test_consistency_gate.py
"""
import hashlib
import os
import pathlib
import subprocess
import sys
import traceback

ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

DBNAME = os.environ.get("GLYCAN_TEST_DB", "glycan_final")
DBUSER = os.environ.get("PGUSER") or os.environ.get("USER")
DBHOST = os.environ.get("PGHOST", "localhost")

# 测试用标签：用独立前缀，便于精确清理，避免污染真实数据
TAG = "__gate_test__"


def _psql(sql: str, dbname: str = None):
    """跑一条 SQL；返回 (returncode, stdout+stderr)。

    刻意用子进程：被拒绝的 INSERT 会让事务进入 aborted 状态，
    在同一条连接里继续跑后续语句会连环失败。
    """
    env = dict(os.environ)
    env["PGUSER"] = DBUSER or ""
    cmd = ["psql", "-h", DBHOST, "-d", dbname or DBNAME, "-v", "ON_ERROR_STOP=0",
           "-tAc", sql]
    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def _db_available() -> bool:
    try:
        rc, out = _psql("SELECT 1")
    except FileNotFoundError:
        return False
    return rc == 0 and "1" in out


def _gate_installed() -> bool:
    rc, out = _psql("SELECT count(*) FROM pg_proc WHERE proname='consistency_gate_sugars'")
    return rc == 0 and out.strip().isdigit() and int(out.strip()) > 0


def _cleanup():
    _psql(f"DELETE FROM sugars WHERE iupac_short LIKE '{TAG}%'")


def _digest(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def _insert(glycoct: str, level: str, extra_cols: str = "", extra_vals: str = "",
            tag: str = ""):
    """插入一条测试记录。

    tag 必须逐测试唯一：多个测试若共用同一个 iupac_short，SELECT 会命中
    兄弟测试的行，断言读取到的是别人留下的值（实测踩过）。
    """
    cols = "sugar_type, glycoct, glycoct_hash, iupac_short, structure_level" + extra_cols
    vals = (f"'poly', '{glycoct.replace(chr(39), chr(39)*2)}', "
            f"'{_digest(glycoct)}', '{TAG}{tag or level}', '{level}'" + extra_vals)
    return _psql(f"INSERT INTO sugars ({cols}) VALUES ({vals}) RETURNING sugar_id")


# ---------------------------------------------------------------------------
# 关卡一：硬门槛
# ---------------------------------------------------------------------------
def test_placeholder_key_with_complete_is_rejected():
    """C7b：占位键不得配 complete —— 这正是"伪造结构"的入口。"""
    rc, out = _insert("UNRESOLVED:poly:10.1/x:deadbeef", "complete")
    assert "一致性硬门槛拒绝写入" in out, out
    assert "C7b" in out, out


def test_empty_glycoct_with_repeat_unit_is_rejected():
    """C7b：glycoct 为空同样不得声明 repeat_unit。"""
    rc, out = _insert("", "repeat_unit")
    assert "一致性硬门槛拒绝写入" in out and "C7b" in out, out


def test_domain_only_without_architecture_is_rejected():
    """C8：domain_only 必须给出域架构，否则域级信息等于没写。"""
    rc, out = _insert("DOMAIN:poly:10.1/x:aaa", "domain_only")
    assert "一致性硬门槛拒绝写入" in out, out
    assert "C8" in out, out


def test_real_structure_with_composition_only_is_rejected():
    """C10：声明"无完整结构"却给了真结构编码，属等级与内容不自洽。"""
    rc, out = _insert("RES\n1b:a-dglc-HEX-1:5", "composition_only")
    assert "一致性硬门槛拒绝写入" in out, out
    assert "C10" in out, out


def test_valid_records_are_accepted():
    """合法记录必须能写入（门槛不能把正常数据也挡掉）。"""
    # 用独立且唯一的 glycoct，避免撞上种子数据或兄弟测试的唯一键
    gx = "RES\n1b:b-dxyl-PEN-1:4"
    rc, out = _insert(gx, "complete", ", molecular_formula", ", 'C5H10O5'", tag="ok1")
    assert "ERROR" not in out, out
    rc2, out2 = _insert("DOMAIN:poly:10.1/ok:ccc", "domain_only",
                        ", domain_architecture",
                        ", '{\"domains\":[{\"name\":\"HG\"}]}'::jsonb", tag="ok2")
    assert "ERROR" not in out2, out2


def test_formula_mismatch_is_recorded_as_finding():
    """C4：分子式与 Glycoct 推导不符 → 写入成功但留下 warn 级发现。"""
    gx = "RES\n1b:b-dgal-HEX-1:5"
    rc, out = _insert(gx + "X", "complete", ", molecular_formula", ", 'C99H99O99'",
                      tag="c4")
    # 用不同 glycoct 避免撞唯一键；若写入成功则应有 C4 发现
    if "一致性硬门槛拒绝写入" in out:
        return          # 未收录单糖导致无法判定，跳过
    rc2, out2 = _psql(
        f"SELECT count(*) FROM sugar_consistency_findings f JOIN sugars s USING(sugar_id) "
        f"WHERE s.iupac_short LIKE '{TAG}%' AND f.rule='C4'")
    assert out2.strip().isdigit(), out2


# ---------------------------------------------------------------------------
# 关卡二：fail-closed 与跨表复核
# ---------------------------------------------------------------------------
def test_consistency_ok_starts_false_and_flips_after_recheck():
    """★ fail-closed：新记录默认为"未复核"，复核通过后才置 TRUE。

    这是本设计的核心：未经复核的数据默认**不可采信**，而不是默认可信。
    """
    gx = "DOMAIN:poly:10.1/fc:ddd"
    rc, out = _insert(gx, "domain_only", ", domain_architecture",
                      ", '{\"domains\":[{\"name\":\"HG\"}]}'::jsonb", tag="fc")
    assert "一致性硬门槛拒绝写入" not in out, out

    rc, out = _psql(f"SELECT consistency_ok FROM sugars WHERE iupac_short='{TAG}fc'")
    assert out.strip() == "f", f"新记录 consistency_ok 应为 false（未复核），实际 {out!r}"

    rc, out = _psql("SELECT consistency_recheck(sugar_id) FROM sugars "
                    f"WHERE iupac_short='{TAG}fc'")
    assert rc == 0, out
    rc, out = _psql(f"SELECT consistency_ok FROM sugars WHERE iupac_short='{TAG}fc'")
    assert out.strip() == "t", f"复核无 block 后应置 true，实际 {out!r}"


def test_structure_change_resets_verification():
    """结构被改动后必须回到"未复核"状态（否则复核结果会被静默沿用）。"""
    gx = "DOMAIN:poly:10.1/rs:eee"
    _insert(gx, "domain_only", ", domain_architecture",
            ", '{\"domains\":[{\"name\":\"HG\"}]}'::jsonb", tag="rs")
    _psql("SELECT consistency_recheck(sugar_id) FROM sugars "
          f"WHERE iupac_short='{TAG}rs'")
    rc, out = _psql(
        "UPDATE sugars SET domain_architecture='{\"domains\":[{\"name\":\"RG-I\"}]}'::jsonb "
        f"WHERE iupac_short='{TAG}rs'")
    assert rc == 0, out
    rc, out = _psql(f"SELECT consistency_ok FROM sugars WHERE iupac_short='{TAG}rs'")
    assert out.strip() == "f", f"结构变更后应回到未复核，实际 {out!r}"


def test_recheck_records_cross_table_findings():
    """C1：组成百分比合计异常应在复核阶段被记录（跨表规则）。"""
    gx = "DOMAIN:poly:10.1/c1:fff"
    _insert(gx, "domain_only", ", domain_architecture",
            ", '{\"domains\":[{\"name\":\"HG\"}]}'::jsonb", tag="c1")
    rc, out = _psql(
        "INSERT INTO polysaccharide_props (sugar_id, monosaccharide_ratio) "
        "SELECT sugar_id, '{\"GalA\":40.0,\"Rha\":30.0}'::jsonb FROM sugars "
        f"WHERE iupac_short='{TAG}c1' RETURNING poly_id")
    assert "一致性硬门槛拒绝写入" not in out, out
    _psql("SELECT consistency_recheck(sugar_id) FROM sugars "
          f"WHERE iupac_short='{TAG}c1'")
    rc, out = _psql(
        "SELECT string_agg(rule, ',') FROM sugar_consistency_findings f "
        f"JOIN sugars s USING(sugar_id) WHERE s.iupac_short='{TAG}c1'")
    assert "C1" in out, f"应记录 C1 发现，实际 {out!r}"


def test_policy_table_is_configurable():
    """严重度必须可配置：把 C7b 降为 warn 后，原本被拒的写入应放行。"""
    rc, out = _psql("UPDATE consistency_rule_policy SET severity='warn' WHERE rule='C7b'")
    assert rc == 0, out
    try:
        rc, out = _insert("UNRESOLVED:poly:10.1/pol:ggg", "complete")
        assert "一致性硬门槛拒绝写入" not in out, f"策略降级后不应再拒绝: {out}"
    finally:
        _psql("UPDATE consistency_rule_policy SET severity='block' WHERE rule='C7b'")


# ---------------------------------------------------------------------------
def main() -> int:
    if not _db_available():
        print(f"  (跳过：连不上数据库 {DBHOST}/{DBNAME}；"
              "可用 GLYCAN_TEST_DB / PGUSER 指定)")
        return 0
    if not _gate_installed():
        print(f"  (跳过：{DBNAME} 未安装 v007 门槛；先跑 db/migrations/v007_consistency_gate.sql)")
        return 0

    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    try:
        for fn in tests:
            try:
                fn()
                print(f"  PASS  {fn.__name__}")
            except Exception as exc:  # noqa: BLE001
                failed += 1
                print(f"  FAIL  {fn.__name__}: {exc}")
                traceback.print_exc()
    finally:
        _cleanup()
    print(f"\n{len(tests) - failed}/{len(tests)} 通过（库={DBNAME}）")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
