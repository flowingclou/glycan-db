#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""回归测试：针对仓库中**真实出现过**的缺陷，防止修复被无意回退。

覆盖范围（全部不依赖数据库 / 网络）:
  1. 结构编码缺失时的占位符必须逐结构唯一（否则不同多糖互相覆盖）
  2. 占位符必须确定性可复现（保证重复导入幂等）
  3. 有真实 GlycoCT 时原样透传
  4. 多糖摘要命名
  5. 首页标题抽取不被页眉期刊名劫持（曾导致所有多糖都叫 unknown polysaccharide）

用法:
  python3 tests/test_regressions.py
"""
import pathlib
import sys
import traceback
import importlib.util

ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from glycan_etl import core                      # noqa: E402
from glycan_etl import table_parser as tp        # noqa: E402
from glycan_etl import glycoct as gx             # noqa: E402


def _load_batch_etl():
    """加载 pipelines/batch_etl.py（不在 glycan_etl 包内，用显式路径加载）。"""
    spec = importlib.util.spec_from_file_location(
        "batch_etl_under_test", str(ROOT / "pipelines" / "batch_etl.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _poly(iupac, doi, residues):
    return core.GlycanRecord(sugar_type="poly", iupac_short=iupac, doi=doi,
                             residues=residues)


# ---------------------------------------------------------------------------
# 1. 结构覆盖缺陷（最严重：不同结构被合并成一条，谱图错挂）
# ---------------------------------------------------------------------------
def test_placeholder_unique_per_structure():
    a = _poly("hawthorn", "10.1/a", [core.Residue(1, "GalA", anomer="a", parent_carbon=4)])
    b = _poly("hawthorn", "10.1/a", [core.Residue(1, "GalA", anomer="a", parent_carbon=3)])
    c = _poly("hawthorn", "10.1/a", [core.Residue(1, "Gal", anomer="a", parent_carbon=4)])
    codes = {core.resolve_glycoct(x) for x in (a, b, c)}
    assert len(codes) == 3, f"不同结构必须得到不同编码, 实际 {codes}"


def test_placeholder_separates_documents():
    """同一占位名 + 不同文献：不得塌缩为同一条记录。"""
    a = _poly("unknown polysaccharide", "10.1/a", [])
    b = _poly("unknown polysaccharide", "10.1/b", [])
    assert core.resolve_glycoct(a) != core.resolve_glycoct(b)


def test_placeholder_is_deterministic():
    def mk():
        return _poly("x", "10.1/a", [core.Residue(1, "Gal", anomer="b", parent_carbon=6)])
    assert core.resolve_glycoct(mk()) == core.resolve_glycoct(mk()), "重复导入必须命中同一编码"


def test_placeholder_never_empty():
    """空编码会撞 UNIQUE 约束并触发 ON CONFLICT 覆盖 —— 必须永不返回空串。"""
    for rec in (_poly("", "10.1/a", []), _poly("x", "10.1/a", []),
                core.GlycanRecord(sugar_type="mono", iupac_short="Glc")):
        code = core.resolve_glycoct(rec)
        assert code and code.strip(), f"占位编码不得为空: {code!r}"


def test_real_glycoct_passthrough():
    real = "RES 1b:b-lglcp-1:5|2:x"
    rec = core.GlycanRecord(sugar_type="mono", iupac_short="β-D-Glcp", glycoct=real)
    assert core.resolve_glycoct(rec) == real


def test_real_glycoct_survives_whitespace():
    rec = core.GlycanRecord(sugar_type="mono", iupac_short="Glc", glycoct="  RES 1b:x  ")
    assert core.resolve_glycoct(rec) == "RES 1b:x"


# ---------------------------------------------------------------------------
# 2. 多糖摘要命名（避免所有多糖都叫同一个占位串）
# ---------------------------------------------------------------------------
def test_poly_summary_name():
    res = [core.Residue(1, "GalA"), core.Residue(2, "GalA"),
           core.Residue(3, "Gal"), core.Residue(4, "Ara")]
    name = tp._poly_summary_name(res)
    assert name.startswith("poly[") and "GalA x2" in name, name
    assert "Gal" in name and "Ara" in name


def test_poly_summary_name_empty():
    assert "unresolved" in tp._poly_summary_name([])


# ---------------------------------------------------------------------------
# 3. 首页标题抽取（曾把页眉期刊名当成标题，连带多糖名提取失效）
# ---------------------------------------------------------------------------
def _page(words):
    return {"page": 1, "n_pages": 1, "text": "", "words": words, "lines": [],
            "width": 595, "height": 842}


def _w(text, x0, x1, top, size):
    return {"text": text, "x0": x0, "x1": x1, "top": top,
            "bottom": top + size, "size": size}


def test_title_picks_largest_font_body_line():
    words = [
        # 页眉期刊名（字号较小）
        _w("Carbohydrate", 40, 120, 20, 8.0),
        _w("Polymers", 125, 170, 20, 8.0),
        # 正文标题（字号最大）
        _w("Structural", 60, 130, 200, 14.0),
        _w("characterisation", 134, 230, 200, 14.0),
        _w("of", 234, 248, 200, 14.0),
        _w("hawthorn", 252, 320, 200, 14.0),
        _w("polysaccharide", 324, 430, 200, 14.0),
    ]
    meta = tp.extract_meta([_page(words)])
    assert meta["title"] and "hawthorn" in meta["title"], meta
    assert meta["polysaccharide"] == "hawthorn", meta


def test_title_rejects_journal_header():
    """期刊名即使字号最大, 也不能被当成标题。"""
    words = [
        _w("International", 40, 120, 20, 14.0),
        _w("Journal", 124, 180, 20, 14.0),
        _w("of", 184, 196, 20, 14.0),
        _w("Biological", 200, 260, 20, 14.0),
        _w("Macromolecules", 264, 360, 20, 14.0),
        _w("A", 60, 70, 300, 11.0),
        _w("real", 74, 100, 300, 11.0),
        _w("paper", 104, 140, 300, 11.0),
        _w("title", 144, 178, 300, 11.0),
    ]
    assert tp.extract_title([_page(words)]) is None


def test_polysaccharide_name_skips_stopwords():
    words = [_w("Study", 60, 100, 200, 14.0),
             _w("on", 104, 122, 200, 14.0),
             _w("a", 126, 134, 200, 14.0),
             _w("polysaccharide", 138, 240, 200, 14.0)]
    meta = tp.extract_meta([_page(words)])
    assert meta["polysaccharide"] is None, meta


# ---------------------------------------------------------------------------
# 4. 标准 GlycoCT 生成（P0-1）
#     旧版手拼的 "RES 1b:a-lglcp-1:5|2:x" 不是合法 GlycoCT，
#     外部工具（glypy / GlyTouCan）无法解析，AI 平台也就无法做结构比对。
# ---------------------------------------------------------------------------
def _res(seq, name, anomer="unknown", parent=None, reducing=False, branch=0):
    return core.Residue(residue_seq=seq, monosaccharide_name=name, ring_form="p",
                        anomer=anomer, is_reducing_end=reducing,
                        parent_carbon=parent, linkage_branch=branch)


def test_glycoct_mono_is_standard_format():
    out = gx.build_glycoct([_res(1, "Glc", anomer="b", reducing=True)], "mono")
    assert out["level"] == "complete"
    assert out["glycoct"] == "RES\n1b:b-dglc-HEX-1:5", out["glycoct"]


def test_glycoct_oligo_linkage_direction():
    """供体提供异头碳、受体提供羟基氧：写法必须是 <受体>o(<位点>+1)<供体>d。"""
    residues = [_res(1, "Glc", anomer="a", parent=1),          # 非还原端，C1 参与键
                _res(2, "Glc", anomer="unknown", parent=4, reducing=True)]  # 还原端 O4
    out = gx.build_glycoct(residues, "oligo")
    assert out["level"] == "complete"
    assert out["glycoct"].splitlines()[-1] == "1:2o(4+1)1d", out["glycoct"]


def test_glycoct_uronic_acid_and_deoxy():
    ga = gx.build_glycoct([_res(1, "GalA", anomer="a")], "mono")["glycoct"]
    assert "|6:a" in ga, ga              # 糖醛酸：C6 为羧酸
    rha = gx.build_glycoct([_res(1, "Rha", anomer="a")], "mono")["glycoct"]
    assert "lman" in rha and "|6:d" in rha, rha   # 6-脱氧-L-甘露糖


def test_glycoct_nacetyl_uses_substituent_residue():
    out = gx.build_glycoct([_res(1, "GlcNAc", anomer="b")], "mono")["glycoct"]
    assert "2s:n-acetyl" in out and "1:1d(2+1)2n" in out, out


def test_glycoct_composition_only_not_fabricated():
    """连接顺序未知的杂多糖不得伪造 GlycoCT。"""
    residues = [_res(1, "GalA", anomer="a", parent=4),
                _res(2, "Gal", anomer="b", parent=None),
                _res(3, "Ara", anomer="a", parent=5)]
    out = gx.build_glycoct(residues, "poly")
    assert out["level"] == "composition_only"
    assert out["glycoct"] is None, out["glycoct"]
    assert out["composition"], out


def test_glycoct_unknown_monosaccharide_degrades_safely():
    out = gx.build_glycoct([_res(1, "Unobtainium", anomer="a")], "mono")
    assert out["glycoct"] is None
    assert out["level"] == "composition_only"


def test_generated_glycoct_passes_glypy():
    """有 glypy 时做闭环校验（无 glypy 则跳过，属可选依赖）。"""
    samples = [
        gx.build_glycoct([_res(1, "Glc", anomer="b")], "mono")["glycoct"],
        gx.build_glycoct([_res(1, "Glc", anomer="a", parent=1),
                          _res(2, "Glc", parent=4, reducing=True)], "oligo")["glycoct"],
        gx.build_glycoct([_res(1, "GlcNAc", anomer="b")], "mono")["glycoct"],
    ]
    for text in samples:
        ok, msg = gx.validate_glycoct(text)
        if ok is None:
            return          # 未安装 glypy，跳过
        assert ok is True, f"{text!r} 校验失败: {msg}"


def test_core_uses_standard_glycoct():
    """core.parse_block 产出的必须是标准编码，而不是旧版手拼串。"""
    rec = core.parse_block(
        "", "β-D-Glcp\n1H NMR (500 MHz, D2O) δ 4.64 (d, J = 7.9 Hz, 1H, H-1).",
        "10.x", "J", 2024)[0]
    assert rec.glycoct and "RES" in rec.glycoct.splitlines()
    assert rec.structure_level == "complete"
    assert not rec.glycoct.startswith("RES 1b:"), rec.glycoct


# ---------------------------------------------------------------------------
# 5. 首页标题 / 多糖命名（多行标题曾只取一行，命名曾完全取不到）
# ---------------------------------------------------------------------------
def test_title_joins_multiline_and_drops_journal_header():
    words = [
        # 页眉期刊名（同为最大字号，必须剔除）
        _w("Journal", 40, 100, 20, 14.0), _w("of", 104, 116, 20, 14.0),
        _w("Pharmaceutical", 120, 200, 20, 14.0), _w("Analysis", 204, 260, 20, 14.0),
        # 跨 3 行的正文标题
        _w("Structural", 60, 130, 180, 14.0),
        _w("characterization", 134, 230, 180, 14.0),
        _w("of", 234, 248, 180, 14.0),
        _w("a", 252, 260, 180, 14.0),
        _w("novel", 264, 300, 180, 14.0),
        _w("polysaccharide", 304, 400, 180, 14.0),
        _w("from", 60, 96, 200, 14.0),
        _w("Polygonati", 100, 170, 200, 14.0),
        _w("Rhizoma", 174, 230, 200, 14.0),
    ]
    meta = tp.extract_meta([_page(words)])
    assert "Structural" in meta["title"] and "Rhizoma" in meta["title"], meta
    assert "Journal" not in meta["title"], meta


def test_polysaccharide_name_from_naming_sentence():
    """标题只给泛称时，从正文命名句取文献内编号（如 SPR-1）。"""
    page = _page([])
    page["text"] = "A purified polysaccharide was obtained, referred to as SPR-1."
    meta = tp.extract_meta([page])
    assert meta["polysaccharide"] == "SPR-1", meta


def test_empty_parse_yields_no_record():
    """没有任何实质数据时不得产出记录 —— 否则报告会把失败显示成成功。"""
    rec = tp.build_record({"polysaccharide": None}, {"name": "x"}, [], [], None)
    assert rec is None


# ---------------------------------------------------------------------------
# 6. 正文连接式 → 完整结构（作者已用 2D NMR 推断好的结论句）
# ---------------------------------------------------------------------------
def test_normalize_structure_text_slash_vs_arrow():
    """'/' 后接 '[' 或数字才是被误映射的箭头；'α/β' 里的斜杠要保留。"""
    assert tp.normalize_structure_text("/4)-b-D-Galp") == "→4)-b-D-Galp"
    assert "a/b" in tp.normalize_structure_text("→4)-a/b-D-Glcp")


def test_parse_sequence_expression_repeat_block_and_acid():
    """重复块 [X]n 展开；糖醛酸 'GalpA' 的 A 需拼回糖名。"""
    residues = tp.parse_sequence_expression("→[4)-b-D-Galp-(1]3→4)-a-D-GalpA-(1")
    names = [r.monosaccharide_name for r in residues]
    assert names == ["Gal", "Gal", "Gal", "GalA"], names
    # 前三个代表重复单元，parent_carbon 都是 4（被取代位点）
    assert all(r.parent_carbon == 4 for r in residues)
    assert residues[0].anomer == "b" and residues[-1].anomer == "a"


def test_parse_sequence_expression_terminal_residue():
    """链末端残基没有 '(1'（不再连出），靠后接标点收尾，并标为还原端。"""
    residues = tp.parse_sequence_expression("→4)-a-D-Glcp-(1→4)-a-D-Glcp, and the")
    assert len(residues) == 2, residues
    assert residues[-1].is_reducing_end is True
    assert residues[0].is_reducing_end is False


def test_extract_linkage_sequence_from_conclusion_sentence():
    text = ("Hence the connection of the main chain was "
            "→[4)-b-D-Galp-(1]2→4,6)-b-D-Galp-(1→4)-a-D-Glcp, and the connection "
            "of the branch chains were R1: b-D-Galp-(1")
    residues = tp.extract_linkage_sequence(text)
    assert [r.monosaccharide_name for r in residues] == ["Gal", "Gal", "Gal", "Glc"]
    assert residues[2].linkage_branch == 6        # 4,6-双取代为分支点


def test_extract_linkage_sequence_ignores_caption():
    """图注 'The (A) main chain, (B) branched-chain ...' 不得被当成结论句。"""
    assert tp.extract_linkage_sequence(
        "Fig. 3. The (A) main chain, (B) branched-chain, and (C) model of SPR-1.") == []


# ---------------------------------------------------------------------------
# 7. 解析出 0 条时的原因分类（区分"文献类型不适用"与"解析未覆盖"）
# ---------------------------------------------------------------------------
def test_explain_zero_records():
    be = _load_batch_etl()
    f = be.explain_zero_records
    assert "扫描件" in f("")
    assert "非结构表征" in f("Molecular dynamics simulation of beta-glucan.")
    assert "无位移数值表" in f("The NMR spectra were computed for this model.")
    shifts = "chemical shift assignments " + " ".join(f"{3.0 + i * 0.11:.2f}" for i in range(30))
    assert "未解析成功" in f(shifts)


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
