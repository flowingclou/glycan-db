#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""域级结论抽取 + 交叉一致性校验 的回归测试。

覆盖两类真实缺陷（均来自山楂多糖 HP，IJBM 2025, 10.1016/j.ijbiomac.2025.145713）：

A. 域级结论被整条丢弃
   该文正文写有 "HP is mainly composed of a large number of HG domains and a
   small number of RG-I domains with side chains."，但这句话**不含连接式箭头**，
   因此 extract_linkage_sequence() 的触发词抓不到，记录被降级成 composition_only，
   "主链是 HG + 侧链在 RG-I 上"这一层信息完全丢失。
   注意：完整结构图在 Fig. 3E，而该面板是**位图**，文本层中不存在 —— 所以
   域级结论句是纯文本管线能拿到的最高一档结构信息。

B. 论文内部数据互相矛盾却完全静默
   组成表（GalA 85.1% / Man 4.506% / Gal 3.983%，合计仅 93.6%）、
   位移归属表（GalA/Gal/Ara/Rha/GlcA，无 Man）、甲基化分析（多出 Glc）、
   摘要（写 α-D-galactose，而表里实测 β-D-Galp）四者互不相容。

全部测试不依赖数据库与网络；有一项 PDF 端到端测试在找不到样本 PDF 时自动跳过。

用法:
  python3 tests/test_domain_and_consistency.py
  GLYCAN_TEST_PDF=/path/to/paper.pdf python3 tests/test_domain_and_consistency.py
"""
import os
import pathlib
import sys
import traceback

ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from glycan_etl import consistency as cs      # noqa: E402
from glycan_etl import core                   # noqa: E402
from glycan_etl import glycoct as gx          # noqa: E402
from glycan_etl import table_parser as tp     # noqa: E402

# ---------------------------------------------------------------------------
# 论文原文片段（从 PDF 文本层逐字摘录）
# ---------------------------------------------------------------------------
DOMAIN_SENTENCE = (
    "The monosaccharide composition, methylation results, and 1D and 2D NMR "
    "analyses of HP indicate that HP is a complex polysaccharide comprising "
    "various monosaccharides. These results suggest that HP is mainly composed "
    "of a large number of HG domains and a small number of RG-I domains with "
    "side chains. The possible structural elements are shown in Fig. 3E.")
METHYLATION_SENTENCE = (
    "Methylation analysis showed that the sugar residue-linking sequence of HP "
    "was mainly 4-GalA (p) and included various terminal sugar residues, such as "
    "t-GalA(p), 4-Glc (p), 2, 4-GalA (p), t-GlcA (p), and 4, 6-GalA (p) (Fig. 1E).")
ABSTRACT_SENTENCE = (
    "The polymer chain is composed of \u03b1-D-galacturonic acid, with branched "
    "chains of \u03b1-D-galactose-\u03b1-L-rhamnose, \u03b1-D-galactose, and "
    "\u03b1-L-arabinose.")

# 山楂文 Table 1（组成表）实测值：合计仅 93.589%，且 Man 在任何残基表中都不存在
HAWTHORN_PCT = {"GalA": 85.10, "Man": 4.506, "Gal": 3.983}

# 山楂文 Table 2（位移归属表）实测残基：无 Man，有 Ara/Rha/GlcA
HAWTHORN_RESIDUES = [
    ("GalA", "a"), ("GalA", "a"), ("GalA", "a"), ("GalA", "a"),
    ("GalA", "a"), ("GalA", "a"), ("GalA", "b"),
    ("Gal", "b"), ("Gal", "b"), ("Gal", "b"),
    ("Ara", "a"), ("Ara", "a"), ("Ara", "a"),
    ("Rha", "a"), ("Rha", "a"),
    ("GlcA", "b"),
]


def _residues():
    return [core.Residue(residue_seq=i + 1, monosaccharide_name=n, anomer=a)
            for i, (n, a) in enumerate(HAWTHORN_RESIDUES)]


def _hawthorn_rec(**kw):
    rec = core.GlycanRecord(
        sugar_type="poly", iupac_short="HP", doi="10.1016/j.ijbiomac.2025.145713",
        residues=_residues())
    rec.poly_props = core.PolysaccharideProps(
        monosaccharide_ratio=dict(HAWTHORN_PCT),
        # 甲基化分析给出的连接类型（含残基表没有的 Glc）
        branching="4-GalA(p), t-GalA(p), 4-Glc(p), t-GlcA(p), 4,6-GalA(p)")
    for k, v in kw.items():
        setattr(rec, k, v)
    return rec


# ---------------------------------------------------------------------------
# A. 域级结论抽取 → domain_only
# ---------------------------------------------------------------------------
def test_extract_domain_architecture_hawthorn():
    """真实结论句 → HG(large) + RG-I(small, 带侧链)。"""
    arch = gx.extract_domain_architecture(DOMAIN_SENTENCE)
    names = [d["name"] for d in arch["domains"]]
    assert names == ["HG", "RG-I"], arch
    hg = arch["domains"][0]
    rgi = arch["domains"][1]
    assert hg["quantity"] == "large" and hg["side_chains"] is False, hg
    assert rgi["quantity"] == "small" and rgi["side_chains"] is True, rgi
    assert "composed of" in arch["evidence"], arch["evidence"]


def test_extract_domain_architecture_negatives():
    """没有域名词的句子不得误判。"""
    for bad in (
        "The main chain was \u21924)-\u03b1-D-GalpA-(1\u2192",
        "Fig. 3E shows the possible structural elements.",
        "HP was composed of galactose and glucose.",
        "",
    ):
        assert gx.extract_domain_architecture(bad)["domains"] == [], bad


def test_build_glycoct_domain_only_does_not_fabricate():
    """给出域架构时必须返回 domain_only 且 glycoct=None，绝不伪造序列。"""
    arch = gx.extract_domain_architecture(DOMAIN_SENTENCE)
    out = gx.build_glycoct(_residues(), "poly", domain_architecture=arch)
    assert out["level"] == "domain_only", out
    assert out["glycoct"] is None, out["glycoct"]
    assert out["composition"], out
    assert any("域级骨架" in w for w in out["warnings"]), out["warnings"]
    # 保留组成式的同时带上域架构
    assert out["domain_architecture"]["domains"], out


def test_build_glycoct_without_domain_unchanged():
    """不传域架构时行为不变（防止影响既有等级判定）。"""
    out = gx.build_glycoct(_residues(), "poly")
    assert out["level"] == "composition_only", out
    assert out["glycoct"] is None


def test_domain_only_key_is_deterministic_and_unique():
    """domain_only 的去重键必须确定性、非空，且不同域架构互不冲突。"""
    arch = gx.extract_domain_architecture(DOMAIN_SENTENCE)
    a = core.GlycanRecord(sugar_type="poly", iupac_short="HP",
                          doi="10.1/x", structure_level="domain_only",
                          domain_architecture=arch)
    b = core.GlycanRecord(sugar_type="poly", iupac_short="HP",
                          doi="10.1/x", structure_level="domain_only",
                          domain_architecture=arch)
    assert core.resolve_glycoct(a) == core.resolve_glycoct(b), "重复导入必须命中同一键"
    assert core.resolve_glycoct(a).startswith("DOMAIN:"), core.resolve_glycoct(a)
    c = core.GlycanRecord(sugar_type="poly", iupac_short="HP", doi="10.1/x",
                          structure_level="domain_only",
                          domain_architecture={"domains": [
                              {"name": "HG", "quantity": "large", "side_chains": False}]})
    assert core.resolve_glycoct(a) != core.resolve_glycoct(c), "不同域架构不得塌缩"
    d = core.GlycanRecord(sugar_type="poly", iupac_short="HP", doi="10.1/y",
                          structure_level="domain_only", domain_architecture=arch)
    assert core.resolve_glycoct(a) != core.resolve_glycoct(d), "不同文献不得塌缩"


# ---------------------------------------------------------------------------
# B. 交叉一致性校验
# ---------------------------------------------------------------------------
def test_canon_mono_does_not_eat_sugar_name_letters():
    """剥离构型前缀时不得吃掉糖名尾字母（曾把 GalA 变成 GlA → 全部识别失败）。"""
    for raw, want in (("GalA", "GalA"), ("GlcA", "GlcA"), ("Gal", "Gal"),
                      ("Glc", "Glc"), ("GlcNAc", "GlcNAc"), ("GalNAc", "GalNAc"),
                      ("Galacturonic acid", "GalA"), ("Glucuronic acid", "GlcA"),
                      ("b-D-Galp", "Gal"), ("\u03b1-D-GalpA", "GalA"),
                      ("GE1,4", "GalA"), ("4MeGlcA", "GlcA"), ("Rha1,2", "Rha")):
        assert cs.canon_mono(raw) == want, f"{raw!r} -> {cs.canon_mono(raw)!r} (期望 {want})"
    assert cs.canon_mono("unknown") is None
    assert cs.canon_mono("") is None


def test_composition_from_residues_accepts_dicts():
    """dry-run 出来的 dict 残基也必须能计数（曾只支持 dataclass）。"""
    rows = [{"monosaccharide_name": "GalA"}, {"monosaccharide_name": "GalA"},
            {"monosaccharide_name": "Rha"}]
    assert cs.composition_from_residues(rows) == {"GalA": 2, "Rha": 1}
    assert cs.composition_from_residues(_residues())["GalA"] == 7


def test_methylation_split_keeps_internal_comma():
    """按逗号切分时不得把 "4,6-GalA(p)" 切成 "6-GalA(p)"。"""
    rec = _hawthorn_rec()
    types = cs._methylation_from(rec)
    assert "4,6-GalA(p)" in types, list(types)
    assert "t-GalA(p)" in types and "4-Glc(p)" in types, list(types)


def test_c1_percent_sum_flagged():
    """组成百分比合计偏离 100% → C1 flagged。"""
    f = cs.check_composition_percent_sum(HAWTHORN_PCT)
    assert len(f) == 1 and f[0]["rule"] == "C1", f
    assert f[0]["severity"] == "flagged"
    assert abs(f[0]["fields"]["total"] - 93.589) < 0.01, f[0]["fields"]
    # 正常值不报
    assert cs.check_composition_percent_sum({"GalA": 60.0, "Gal": 40.0}) == []


def test_c2_component_set_mismatch():
    """组成表有 Man、残基表没有；残基表有 Ara/Rha/GlcA、组成表没有。"""
    f = cs.check_component_sets(HAWTHORN_PCT, _residues())
    assert len(f) == 1 and f[0]["rule"] == "C2", f
    assert f[0]["fields"]["only_in_percentages"] == ["Man"], f[0]["fields"]
    assert f[0]["fields"]["only_in_residues"] == ["Ara", "GlcA", "Rha"], f[0]["fields"]


def test_c3_methylation_only_sugar():
    """甲基化分析出现位移表没有的 Glc → C3。"""
    rec = _hawthorn_rec()
    f = cs.check_component_sets(HAWTHORN_PCT, _residues(),
                                methylation=cs._methylation_from(rec))
    rules = {x["rule"] for x in f}
    assert "C3" in rules, f
    c3 = [x for x in f if x["rule"] == "C3"][0]
    assert c3["fields"]["only_in_methylation"] == ["Glc"], c3["fields"]


def test_c5_anomer_conflict_abstract_vs_table():
    """摘要写 α-D-galactose，归属表 Gal 实测全是 β → C5 flagged。"""
    f = cs.check_anomer_conflict(ABSTRACT_SENTENCE, _residues())
    assert len(f) == 1 and f[0]["rule"] == "C5", f
    assert f[0]["fields"]["mono"] == "Gal", f[0]["fields"]
    assert f[0]["fields"]["stated"] == ["a"], f[0]["fields"]
    assert f[0]["fields"]["observed"] == ["b"], f[0]["fields"]
    # 构型一致时不得误报
    ok = [core.Residue(1, "Gal", anomer="a")]
    assert cs.check_anomer_conflict(ABSTRACT_SENTENCE, ok) == []


def test_c4_formula_vs_residues():
    """分子式对账：C12H22O11 ↔ 两个己糖成苷 = 一致；C6H12O6 则不符。"""
    good = cs.check_formula_vs_residues("C12H22O11", "Glc2")
    assert good == [], good
    bad = cs.check_formula_vs_residues("C6H12O6", "Glc2")
    assert len(bad) == 1 and bad[0]["rule"] == "C4", bad
    # N-乙酰氨基糖的 N 必须算进去（曾经用 Glc 的分子式冒充 GlcNAc）
    calc = cs.formula_from_composition("GlcNAc", n_bonds=0)
    assert calc.get("N") == 1, calc
    # 结构不完整时只给 warning，不误判为解析错误
    soft = cs.check_formula_vs_residues("C6H12O6", "GalA7,Ara3,Gal3,Rha2,GlcA",
                                        structure_level="domain_only")
    assert soft and soft[0]["severity"] == "warning", soft


def test_c6_methylation_vs_residue_counts():
    rec = _hawthorn_rec()
    f = cs.check_methylation_counts(cs._methylation_from(rec), _residues())
    assert len(f) == 1 and f[0]["rule"] == "C6", f


def test_cross_check_end_to_end_on_hawthorn():
    """一次跑完：C1/C2/C3/C5 全部命中；apply_to_record 置 flagged 并写 notes。"""
    rec = _hawthorn_rec()
    out = cs.cross_check(rec, statement_text=ABSTRACT_SENTENCE)
    rules = {f["rule"] for f in out["findings"]}
    assert {"C1", "C2", "C3", "C5"} <= rules, rules
    assert out["flagged"] is True and out["n_flagged"] >= 4, out["summary"]
    n = cs.apply_to_record(rec, out["findings"])
    assert n == len(out["findings"]), (n, len(out["findings"]))
    assert rec.qc_status == "flagged", rec.qc_status
    assert any("一致性[C1]" in x for x in rec.qc_notes), rec.qc_notes


def test_cross_check_clean_record_not_flagged():
    """干净记录不得被误报（避免把无冲突数据也标 flagged）。"""
    rec = core.GlycanRecord(
        sugar_type="poly", iupac_short="clean", doi="10.1/z",
        residues=[core.Residue(1, "Glc", anomer="a"),
                  core.Residue(2, "Glc", anomer="a")],
        structure_level="repeat_unit", composition="Glc2")
    rec.poly_props = core.PolysaccharideProps(
        monosaccharide_ratio={"Glc": 100.0}, branching="4-Glc(p)")
    out = cs.cross_check(rec)
    assert out["findings"] == [], out["summary"]


def test_consistency_flag_survives_build_record():
    """build_record 不得把一致性校验给出的 flagged 覆盖成 dry-run-ok/passed。

    实测踩过：qc_status 的基础赋值写在一致性校验**之后**，导致报告里
    flagged 恒为 0（校验明明命中了 C1/C2/C3/C5）。
    """
    meta = {"doi": "10.1/hp", "journal": "J", "year": 2025,
            "polysaccharide": "HP", "title": "t"}
    kv = {"name": "HP", "monosaccharide_mol_pct": dict(HAWTHORN_PCT)}
    entries = [{"code": f"R{i}", "iupac": "\u21924)-\u03b1-D-GalpA-(1\u2192",
                "positions": {"1": {"H": 5.0, "C": 100.0}, "4": {"H": 4.3, "C": 78.0}}}
               for i in range(3)]
    rec = tp.build_record(meta, kv, entries, [], 7, full_text="")
    assert rec is not None, "应产出记录"
    assert rec.consistency_findings, "一致性校验应产出发现"
    assert rec.qc_status == "flagged", (
        f"一致性命中后 qc_status 应为 flagged，实际 {rec.qc_status!r}"
        "（基础状态赋值覆盖了校验结果）")
    assert any("一致性[C" in n for n in rec.qc_notes), rec.qc_notes


# ---------------------------------------------------------------------------
# C. 端到端（需要样本 PDF；找不到则跳过）
# ---------------------------------------------------------------------------
def _find_pdf():
    env = os.environ.get("GLYCAN_TEST_PDF")
    cands = [env] if env else []
    cands += [
        "/tmp/glycan_pdf/hawthorn_main.pdf",
        str(ROOT / "tests" / "fixtures" / "hawthorn_main.pdf"),
    ]
    for c in cands:
        if c and os.path.isfile(c):
            return c
    return None


def test_pdf_end_to_end_domain_only():
    pdf = _find_pdf()
    if not pdf:
        print("      (跳过：未找到样本 PDF，可用 GLYCAN_TEST_PDF=... 指定)")
        return
    d = tp.run(pdf, verbose=False)
    assert d["records"], "应解析出记录"
    rec = d["records"][0]
    assert rec.get("structure_level") == "domain_only", rec.get("structure_level")
    assert rec.get("glycoct") is None, "domain_only 不得生成 GlycoCT"
    names = [x["name"] for x in (rec.get("domain_architecture") or {}).get("domains", [])]
    assert names == ["HG", "RG-I"], names
    notes = " ".join(rec.get("qc_notes") or [])
    for rule in ("C1", "C2", "C3", "C5"):
        assert f"一致性[{rule}]" in notes, (rule, notes[:300])
    assert rec.get("composition", "").startswith("GalA7"), rec.get("composition")


# ---------------------------------------------------------------------------
def main() -> int:
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in tests:
        try:
            fn()
            print(f"  PASS  {fn.__name__}")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"  FAIL  {fn.__name__}: {exc}")
            traceback.print_exc()
    print(f"\n{len(tests) - failed}/{len(tests)} 通过")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
