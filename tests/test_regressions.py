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

ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from glycan_etl import core                      # noqa: E402
from glycan_etl import table_parser as tp        # noqa: E402


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
