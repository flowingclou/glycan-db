#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
batch_etl.py — glycan_etl_v3 批量执行器 v1.0（2026-09-12）

把「单 PDF → 手动改 config → dry-run → 入库」的多步流程，
封装为一次命令，支持整目录批量执行 + 自动汇总报告。

用法:
  ① 批量 dry-run（推荐先用，不连库、零风险）:
       python batch_etl.py --dir /path/to/pdf_folder [--filter SI]
       python batch_etl.py --pdf a.pdf b.pdf c.pdf
     → 逐个解析，汇总报告落在  output/etl_reports/ 下（CSV + Markdown）

  ② 批量入库（PostgreSQL）:
       python batch_etl.py --dir /path/to/pdf_folder --config etl_config.yaml
       python batch_etl.py --pdf a.pdf --config cfg.yaml [--meta meta.yaml]
     → --meta 按文件名映射每份的 doi/journal/year（批量多文献必备），
       缺省则所有 PDF 用 config 里的统一 doi/journal/year

  ③ 自检:
       python batch_etl.py --self-test

依赖: pdfplumber psycopg2-binary pyyaml（与 glycan_etl_v3.py 相同）
"""
import argparse
import csv
import datetime
import importlib.util
import json
import logging
import os
import re
import sys

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("batch_etl")

REPORT_DIR = "etl_reports"

# 将仓库根加入 sys.path，使 `from glycan_etl import core` 在整个仓库中可用
# （本文件位于仓库根下的 pipelines/ 子目录）
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)


# ---------------------------------------------------------------------------
# 加载 glycan_etl_v3 主脚本（同目录优先，可用 --etl 显式指定路径）
# ---------------------------------------------------------------------------
def load_etl(etl_path: str = None):
    """加载解析引擎：优先 glycan_etl 包（仓库模块化方式）。

    也兼容旧式用法：--etl 显式指定散装 glycan_etl_v3.py 路径。
    """
    if etl_path is None:
        try:
            from glycan_etl import core
            return core
        except Exception as exc:  # noqa: BLE001
            log.warning("import glycan_etl.core 失败（%s），回退同目录旧脚本", exc)
        here = os.path.dirname(os.path.abspath(__file__))
        cand = [os.path.join(here, "glycan_etl_v3.py"),
                os.path.join(here, "glycan_etl_v2.py"),
                os.path.join(here, "glycan_etl.py")]
        p = next((c for c in cand if os.path.exists(c)), None)
        if p is None:
            sys.exit("找不到 glycan_etl 包或同目录 glycan_etl*.py，请用 --etl 指定路径")
    elif os.path.exists(etl_path):
        p = etl_path
    else:
        sys.exit(f"--etl 指定的脚本不存在: {etl_path}")
    spec = importlib.util.spec_from_file_location("glycan_etl_mod", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# PDF 收集
# ---------------------------------------------------------------------------
def collect_pdfs(pdf_files=None, pdf_dir=None, filter_kw=None):
    paths = []
    if pdf_files:
        for f in pdf_files:
            if not os.path.exists(f):
                log.warning("跳过不存在的文件: %s", f)
                continue
            paths.append(f)
    if pdf_dir:
        if not os.path.isdir(pdf_dir):
            sys.exit(f"目录不存在: {pdf_dir}")
        for root, _, files in os.walk(pdf_dir):
            for fn in sorted(files):
                if fn.lower().endswith(".pdf"):
                    p = os.path.join(root, fn)
                    if filter_kw and filter_kw.lower() not in os.path.basename(p).lower():
                        continue
                    paths.append(p)
    # 去重保序
    seen, uniq = set(), []
    for p in paths:
        ap = os.path.abspath(p)
        if ap not in seen:
            seen.add(ap)
            uniq.append(p)
    return uniq


# ---------------------------------------------------------------------------
# 单份 PDF 解析（复用主脚本的解析管线, 不连库）
# ---------------------------------------------------------------------------
def parse_one_pdf(gle, pdf_path, doi, journal, year):
    text = gle.extract_pdf_text(pdf_path)
    blocks = gle.split_blocks(text)
    records = []
    for header, body in blocks:
        recs = gle.parse_block(header, body, doi, journal, year)
        records.extend(recs)
    return text, blocks, records


# ===========================================================================
# 自动分流解析层 (auto-dispatch v1.0)
#  对每份 PDF 先轻量探测类型（表格型 / 行内型 / 混合型），再自动选路：
#    - 表格型  -> glycan_etl/table_parser.py（位移归属表 / 分子量表）
#    - 行内型  -> glycan_etl/core.py（正文分子式 NMR 写法）
#    - 混合型  -> 两路都跑并合并去重
#  主路结果为空时自动补跑另一路（fallback），保证不漏。
# ===========================================================================

# 探测规则配置（可调）：
#   判定优先级：table_score>=2 且 table>=inline -> table；
#               inline>=2 且 inline>table      -> inline；
#               两者都命中                     -> mixed；
#               都未命中                       -> inline（默认 core，零风险）
DETECT_CFG = {
    "head_pages": 5,                 # 读取前几页做轻量探测
    # 前几页无任何清晰信号时的扩展探测页数：多糖文献的结构表征数据
    # （甲基化表、位移归属表）常排在第 6~10 页，只读前 5 页会误判为
    # 行内型，进而触发"先跑行内、失败再跑表格"的双跑回退。
    "deep_head_pages": 12,
    # 解析出 0 条时，为判定原因而读取的页数（要覆盖全文才看得准）
    "deep_zero_pages": 30,
    # 表格型特征
    "table_caption_re": r"(?m)^\s*Table\s+\d+[\.:\s]",
    "ppm_token_re": r"\d+\.\d+",
    "col_ppm_min": 3,                # 一行至少这么多个 ppm 数值算"表格数据行"
    "ppm_value_range": (0.3, 230.0),  # 合理化学位移数值范围
    "usec_head_hint": ("1H", "13C", "δ", "ppm"),   # 归属表头提示词
    # 行内型特征（正文分子式 NMR 写法）
    "inline_j_re": r"\bδ\s*\d+\.\d+\s*\([a-z]+\s*,?\s*(J\s*=\s*\d+(\.\d+)?\s*Hz)?",
    "inline_shift_re": r"\d+\.\d{2}\s*\((s|d|t|q|m|br\s*s|br\s*d)\b[^)]*\)",
    "inline_atom_re": r"\bH-\d+|\bC-\d+",
}


def peek_pdf_text(pdf_path: str, head_pages: int = None):
    """读取 PDF 前几页纯文本（含标题/引言/表格区域），用于类型探测。"""
    head_pages = head_pages if head_pages is not None else DETECT_CFG["head_pages"]
    try:
        import pdfplumber
    except ImportError:
        log.warning("缺少 pdfplumber，无法探测 PDF 类型，将按默认行内型解析")
        return None
    parts = []
    try:
        with pdfplumber.open(pdf_path) as pdf:
            for page in pdf.pages[:head_pages]:
                parts.append(page.extract_text() or "")
    except Exception as e:  # noqa: BLE001
        log.warning("探测首页文本失败（%s），将按默认行内型解析: %s", os.path.basename(pdf_path), e)
        return None
    if not parts:
        return None
    return "\n".join(parts)


def _score_table_like(text: str) -> dict:
    """启发式：是否存在化学位移归属表 / 分子量表特征。"""
    if not text:
        return {"captions": 0, "columns": 0, "head_hits": 0, "table_hits": 0}
    cfg = DETECT_CFG
    captions = len(re.findall(cfg["table_caption_re"], text))
    ppm_re = re.compile(r"\d+\.\d+")
    lo, hi = cfg["ppm_value_range"]
    col_rows = 0          # 单行 ≥3 个数值 → 疑似表格数据行
    table_rows = 0        # 数值个数在合理范围内的表格行（排除正文长句）
    for line in text.splitlines():
        nums = [float(x) for x in ppm_re.findall(line)]
        if len(nums) < cfg["col_ppm_min"]:
            continue
        in_range = [n for n in nums if lo <= n <= hi]
        if len(in_range) >= cfg["col_ppm_min"]:
            col_rows += 1
        if len(in_range) >= 2:
            table_rows += 1
    ll = text.lower()
    head_hits = sum(1 for k in cfg["usec_head_hint"] if k.lower() in ll)
    # 判定：多处 Table 提注 + 列式数值行/归属表头提示
    table_hits = 0
    if captions >= 1 and (col_rows >= 2 or head_hits >= 2):
        table_hits += 1
    if captions >= 1 and table_rows >= 2 and head_hits >= 1:
        table_hits += 1
    if col_rows >= 3 and head_hits >= 2:
        table_hits += 1
    return {"captions": captions, "columns": col_rows, "head_hits": head_hits,
            "table_hits": table_hits}


def _score_inline_like(text: str) -> dict:
    """启发式：正文分子式 NMR 写法（如 δ 5.35 (d, J = 7.8) / 3.95 (s, 3H, H-1)）。"""
    if not text:
        return {"shifts": 0, "j": 0, "atoms": 0, "inline_hits": 0}
    cfg = DETECT_CFG
    j = len(re.findall(cfg["inline_j_re"], text))
    shifts = len(re.findall(cfg["inline_shift_re"], text))
    atoms = len(re.findall(cfg["inline_atom_re"], text))
    inline_hits = 0
    if j >= 1:
        inline_hits += 1
    if shifts >= 1:
        inline_hits += 1
    if atoms >= 2:
        inline_hits += 1
    return {"shifts": shifts, "j": j, "atoms": atoms, "inline_hits": inline_hits}


def detect_pdf_type(text: str = None, pdf_path: str = None):
    """轻量探测 PDF 类型（表格型/行内型/混合型）。

    优先级: table >= inline -> table；inline > table -> inline；
    两者都命中 -> mixed；都未命中 -> inline（fallback，默认 core）。
    返回 (mode, evidence)。
    """
    if text is None and pdf_path is not None:
        text = peek_pdf_text(pdf_path)
    if not text:
        return "inline", {"mode": "inline", "fallback": "no-text",
                          "table_hits": 0, "inline_hits": 0}
    ts = _score_table_like(text)
    ns = _score_inline_like(text)
    t, i = ts["table_hits"], ns["inline_hits"]
    evidence = {**ts, **ns}
    if t >= 2 and t >= i:
        mode = "table"
    elif i >= 2 and i > t:
        mode = "inline"
    elif t >= 1 or i >= 1:
        mode = "mixed"
    else:
        mode = "inline"
        evidence["fallback"] = "no-clear-signal"
    evidence["mode"] = mode
    return mode, evidence


# ---------------------------------------------------------------------------
# 表格型记录 -> core.GlycanRecord 归一化（对齐入库前格式）
# ---------------------------------------------------------------------------
def _dataclass_from(blob, cls):
    """按 dataclass 字段白名单从 dict 构造，忽略多余键（兼容不同宽度输出）。"""
    fields = set(cls.__dataclass_fields__)
    return cls(**{k: v for k, v in blob.items() if k in fields})


def table_blob_to_records(blob: dict, doi, journal, year):
    """把 table_parser 输出的 asdict 记录结构归一化为 core.GlycanRecord。"""
    from glycan_etl import core as _core
    rec = _core.GlycanRecord(
        sugar_type=blob.get("sugar_type") or "poly",
        iupac_short=blob.get("iupac_short"),
        glycoct=blob.get("glycoct"),
        molecular_formula=blob.get("molecular_formula"),
        molecular_weight=blob.get("molecular_weight"),
        anomer=blob.get("anomer"),
        doi=doi if doi not in (None, "UNKNOWN") else blob.get("doi"),
        journal=journal if journal is not None else blob.get("journal"),
        year=year if year is not None else blob.get("year"),
        nmr_page=blob.get("nmr_page"),
        structure_level=blob.get("structure_level"),
        composition=blob.get("composition"),
    )
    rec.residues = [_dataclass_from(r, _core.Residue) for r in blob.get("residues") or []]
    if blob.get("physicochemical"):
        rec.physicochemical = _dataclass_from(blob["physicochemical"], _core.Physicochemical)
    if blob.get("poly_props"):
        rec.poly_props = _dataclass_from(blob["poly_props"], _core.PolysaccharideProps)
    for e in blob.get("experiments") or []:
        exp = _core.NMRExperiment(
            nucleus=e.get("nucleus", "1H"),
            experiment_2d=e.get("experiment_2d"),
            solvent=e.get("solvent") or "D2O",
            frequency=e.get("frequency"),
            temperature=e.get("temperature"),
            ph=e.get("ph"),
        )
        exp.peaks = [_dataclass_from(p, _core.Peak1D) for p in e.get("peaks") or []]
        exp.peaks_2d = [_dataclass_from(p, _core.Peak2D) for p in e.get("peaks_2d") or []]
        rec.experiments.append(exp)
    rec.qc_notes = list(blob.get("qc_notes") or [])
    # 归一化 QC 状态：table-parser 的 dry-run-ok 视为通过（已有结构/位移证据）
    qs = blob.get("qc_status", "dry-run-ok")
    rec.qc_status = {"dry-run-ok": "passed"}.get(qs, qs)
    if rec.qc_status == "passed" and not any("table-parser" in n for n in rec.qc_notes):
        rec.qc_notes.append("auto-dispatch: table-parser")
    return rec


def table_dryrun_to_records(d: dict, doi, journal, year):
    """解析 table_parser.run() 返回 dict，产出 core 记录列表。"""
    recs = []
    for blob in d.get("records") or []:
        recs.append(table_blob_to_records(blob, doi, journal, year))
    return recs


def merge_records(records):
    """按 (sugar_type, iupac_short 规范化) 去重合并。

    同一结构（如 16 残基 poly）两路都可能命中，保留信息更全的一条：
    table-parser 记录优先（含残基/位移归属表）；同名行内记录去重。
    """
    key2i = {}
    order = []
    for idx, r in enumerate(records):
        key = (r.sugar_type or "?", (r.iupac_short or "").strip().lower())
        # 同 key 已存在：table-parser 优先；已是最优则保留先到者
        if key in key2i:
            keep = key2i[key]
            cur_better = "table-parser" in (r.qc_notes or [])
            old_better = "table-parser" in (records[keep].qc_notes or [])
            if cur_better and not old_better:
                key2i[key] = idx
            continue
        key2i[key] = idx
        order.append(idx)
    return [records[i] for i in sorted(order)]


def dispatch_pdf(gle, pdf_path, doi, journal, year):
    """自动分流解析：探测类型 -> 选路（/合并）-> 统一输出 core 记录。

    返回 (mode, evidence, text, blocks, records)。
    text/blocks 仅在行内（core）路径非空，表格路径为空（不影响入库与报告）。
    """
    # 旧式散装引擎（非 glycan_etl 包）无数据类，回退原始 core-only 路径
    if not getattr(gle, "GlycanRecord", None):
        text, blocks, records = parse_one_pdf(gle, pdf_path, doi, journal, year)
        return "inline", {"mode": "inline", "fallback": "legacy-etl"}, text, blocks, records

    probe_text = peek_pdf_text(pdf_path)
    mode, ev = detect_pdf_type(text=None, pdf_path=pdf_path) if probe_text is None \
        else detect_pdf_type(text=probe_text)
    # 自适应深探测：前几页没有任何清晰信号时，扩展到更多页重探。
    # 不这样做会把"数据排在 6~10 页"的文献误判为行内型，触发双跑回退
    # （既慢，报告里的 mode 也会误导使用者）。
    if not ev.get("table_hits") and not ev.get("inline_hits"):
        deep = peek_pdf_text(pdf_path, head_pages=DETECT_CFG["deep_head_pages"])
        if deep:
            mode2, ev2 = detect_pdf_type(text=deep)
            if ev2.get("table_hits") or ev2.get("inline_hits"):
                mode, ev = mode2, ev2
                ev["probe"] = "deep"
    text = blocks = None
    records = []

    def _run_inline():
        t = gle.extract_pdf_text(pdf_path)
        bl = gle.split_blocks(t)
        rr = []
        for header, body in bl:
            rr.extend(gle.parse_block(header, body, doi, journal, year))
        return t, bl, rr

    def _run_table():
        try:
            from glycan_etl import table_parser as tp
            d = tp.run(pdf_path, verbose=False)
            return d.get("records"), None, table_dryrun_to_records(d, doi, journal, year)
        except SystemExit:
            raise
        except Exception as e:  # noqa: BLE001
            log.warning("table-parser 解析失败（%s）: %s", os.path.basename(pdf_path), e)
            return [], None, []

    if mode == "table":
        _recs_raw, _, trecs = _run_table()
        if trecs:
            records = trecs
        else:  # fallback：表格主路 0 条 -> 补跑行内 core 合并
            t, bl, irecs = _run_inline()
            text, blocks = t, bl
            records = merge_records(irecs + trecs)
            ev["fallback"] = "table-empty-run-inline"
    elif mode == "inline":
        t, bl, irecs = _run_inline()
        text, blocks = t, bl
        records = irecs
        if not irecs:  # fallback：行内主路 0 条 -> 补跑表格路合并
            _r, _, trecs = _run_table()
            records = merge_records(irecs + trecs)
            ev["fallback"] = "inline-empty-run-table"
    else:  # mixed：优先表格路；仅当确有行内特征时才补跑行内路
        _r, _, trecs = _run_table()
        if trecs and not ev.get("inline_hits"):
            # 行内特征为 0（如整篇 NMR 数据只出现在表格里）：再跑一次全文行内
            # 解析纯属浪费——这条路对 15 页 PDF 要额外数秒，且必然 0 结果。
            records = trecs
            ev["note"] = "table-only (inline 特征为 0)"
        else:
            t, bl, irecs = _run_inline()
            text, blocks = t, bl
            records = merge_records(irecs + trecs)
    return mode, ev, text, blocks, records


# ---------------------------------------------------------------------------
# 报告生成
# ---------------------------------------------------------------------------
def explain_zero_records(text: str) -> str:
    """解析出 0 条时给出原因分类。

    区分"这份文献本来就不含结构表征数据"（分子模拟/机理/综述类）与
    "含数据但当前解析未覆盖"—— 前者无需处理，后者才值得投入改进。
    旧报告对两者都只说"请检查文本层质量或结构写法"，把用户引向错误方向。
    """
    if not text:
        return "无文本层（可能为扫描件，需先 OCR）"
    shift_vals = re.findall(r"\b\d{1,3}\.\d{1,2}\b", text)
    has_shift_head = bool(re.search(r"chemical\s*shift|δ\s*\d|ppm", text, re.I))
    has_nmr = bool(re.search(r"\bNMR\b", text))
    if has_shift_head and len(shift_vals) >= 20:
        return ("含位移/化学位移数据但未解析成功（版式或写法未覆盖，"
                "建议反馈该样本）")
    if has_nmr:
        return "提及 NMR 但无位移数值表（多为分子模拟/机理/综述类，非结构表征）"
    if len(shift_vals) < 10:
        return "无化学位移数据（非结构表征类文献，如分子模拟/机理研究）"
    return "有数值但未识别为结构表征（版式未覆盖，建议反馈该样本）"


def dump_report(items, out_dir):
    """items: list of dict（每份 PDF 一行的汇总）。返回 (csv_path, md_path)。"""
    os.makedirs(out_dir, exist_ok=True)
    st = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_path = os.path.join(out_dir, f"manifest_{st}.csv")
    md_path = os.path.join(out_dir, f"manifest_{st}.md")
    fields = ["file", "mode", "records", "flagged", "passed",
              "types", "qc_issues", "detail"]
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for it in items:
            w.writerow(it)
    # Markdown 摘要
    ok = sum(1 for it in items if it["records"] > 0)
    total = sum(it["records"] for it in items)
    modes = ",".join(sorted({str(it.get("mode", "-")) for it in items}))
    lines = [
        "# glycan_etl 批量执行报告",
        f"生成时间: {datetime.datetime.now().isoformat(timespec='seconds')}",
        "",
        f"- 处理 PDF 份数: **{len(items)}**",
        f"- 识别出记录的总 PDF 份数: **{ok}** / {len(items)}",
        f"- 共识别结构-谱图对: **{total}** 条",
        f"- 自动分流模式: **{modes or '-'}**（table=表格型 / inline=行内型 / mixed=双路合并）",
        "",
        "| 文件 | 分流模式 | 记录数 | flagged | passed | 类型分布 | QC问题 | 明细 |",
        "|------|---------|-------:|-------:|-------:|---------|--------|------|",
    ]
    for it in items:
        lines.append(
            f"| {os.path.basename(it['file'])} | {it.get('mode', '-')} | {it['records']} | "
            f"{it['flagged']} | {it['passed']} | {it['types'] or '-'} | "
            f"{it['qc_issues'] or '-'} | {(it['detail'] or '-')[:120]} |"
        )
    lines.append("")
    lines.append("> 0 条的 PDF 请先看「明细」列的原因分类：")
    lines.append("> - 标注「非结构表征类文献」→ 该文献本就不含结构-谱图数据，无需处理；")
    lines.append("> - 标注「未解析成功 / 版式未覆盖」→ 才是需要改进的样本，建议反馈。")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    return csv_path, md_path


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="glycan_etl 批量执行器")
    ap.add_argument("--pdf", nargs="+", help="一个或多个 PDF 路径")
    ap.add_argument("--dir", help="扫描目录（递归找 *.pdf）")
    ap.add_argument("--filter", default=None, help="仅处理文件名含该关键词的 PDF（配合 --dir）")
    ap.add_argument("--config", default="etl_config.yaml", help="入库配置（含 db 连接）")
    ap.add_argument("--meta", default=None, help="per-file 元数据 yaml: 文件名→doi/journal/year")
    ap.add_argument("--etl", default=None, help="glycan_etl_v3.py 路径（默认同目录自动找）")
    ap.add_argument("--report-dir", default=REPORT_DIR, help="报告输出目录（默认 ./etl_reports）")
    ap.add_argument("--dry-run", action="store_true", help="只解析不连库（默认模式）")
    ap.add_argument("--save-json", action="store_true", help="每份 PDF 额外保存完整解析 JSON")
    ap.add_argument("--self-test", action="store_true", help="内置样例自检（不连库）")
    ap.add_argument("--skip-flagged", action="store_true",
                    help="质控违规记录不入库（默认入库并标记 flagged，便于回溯）")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return
    if not args.pdf and not args.dir:
        ap.error("需要 --pdf 或 --dir")
    if not args.dry_run and not os.path.exists(args.config):
        log.warning("未找到 --config 文件 %s，将按 dry-run 模式执行", args.config)
        args.dry_run = True

    gle = load_etl(args.etl)

    # 配置 / per-file 元数据
    config = {}
    if os.path.exists(args.config):
        import yaml
        with open(args.config, encoding="utf-8") as f:
            config = yaml.safe_load(f) or {}
    default_doi = config.get("doi", "UNKNOWN")
    default_journal = config.get("journal")
    default_year = config.get("year")
    meta = {}
    if args.meta and os.path.exists(args.meta):
        import yaml
        with open(args.meta, encoding="utf-8") as f:
            meta = yaml.safe_load(f) or {}

    pdfs = collect_pdfs(args.pdf, args.dir, args.filter)
    if not pdfs:
        log.warning("未收集到任何 PDF（检查 --dir / --filter / --pdf）")
        return
    log.info("共 %d 份 PDF 待处理", len(pdfs))

    items = []
    out_dir = os.path.abspath(args.report_dir)
    # 提前建目录：--save-json 在本循环内写文件，晚于 dump_report() 的 makedirs
    os.makedirs(out_dir, exist_ok=True)
    for i, pdf in enumerate(pdfs, 1):
        base = os.path.basename(pdf)
        m = meta.get(base) or {}
        doi = m.get("doi", default_doi)
        journal = m.get("journal", default_journal)
        year = m.get("year", default_year)
        log.info("[%d/%d] %s (doi=%s)", i, len(pdfs), base, doi)
        try:
            mode, evidence, text, blocks, records = dispatch_pdf(gle, pdf, doi, journal, year)
        except Exception as e:  # noqa: BLE001
            log.error("解析失败: %s (%s)", base, e)
            items.append({"file": pdf, "mode": "-", "records": 0, "flagged": 0, "passed": 0,
                          "types": "-", "qc_issues": f"ERROR: {e}", "detail": ""})
            continue
        fb = evidence.get("fallback")
        if fb:
            log.info("  [%s] 分流 %s（fallback: %s）", base, mode, fb)
        else:
            log.info("  [%s] 分流 %s", base, mode)
        flagged = sum(1 for r in records if r.qc_status == "flagged")
        passed = sum(1 for r in records if r.qc_status == "passed")
        types = ",".join(sorted({r.sugar_type for r in records}))
        issues = "; ".join(" | ".join(r.qc_notes) for r in records if r.qc_notes)
        detail = ", ".join(f"{r.iupac_short}({r.sugar_type},{len(r.experiments)}exp)"
                           for r in records)
        if not records:
            # 0 条时说明原因：是"文献类型不适用"还是"解析未覆盖"
            detail = explain_zero_records(
                peek_pdf_text(pdf, head_pages=DETECT_CFG["deep_zero_pages"]) or "")
        items.append({
            "file": pdf, "mode": mode, "records": len(records), "flagged": flagged,
            "passed": passed, "types": types, "qc_issues": issues[:400], "detail": detail,
        })
        if args.save_json:
            with open(os.path.join(out_dir, f"{os.path.splitext(base)[0]}_records.json"),
                      "w", encoding="utf-8") as f:
                json.dump([{"iupac_short": r.iupac_short, "sugar_type": r.sugar_type,
                            "anomer": r.anomer, "n_experiments": len(r.experiments),
                            "n_2d_peaks": sum(len(e.peaks_2d) for e in r.experiments),
                            "structure_level": r.structure_level,
                            "composition": r.composition,
                            "glycoct": r.glycoct,
                            "residues": [asdict_min(r) for r in r.residues],
                            "qc_status": r.qc_status, "qc_notes": r.qc_notes}
                           for r in records], f, ensure_ascii=False, indent=2)

        # 入库模式
        if not args.dry_run:
            try:
                import psycopg2
                db = config.get("db")
                if not db:
                    raise RuntimeError("config 中缺少 db 连接信息")
                conn = psycopg2.connect(**db)
                n, skipped = 0, 0
                for rec in records:
                    if rec.qc_status == "flagged" and args.skip_flagged:
                        log.warning("  跳过 flagged 记录: %s (%s)",
                                    rec.iupac_short, rec.qc_notes)
                        skipped += 1
                        continue
                    if rec.qc_status == "flagged":
                        # 默认保留可疑数据并打 flagged 标记：真值库里丢弃比保留
                        # 风险更高（无法回溯），由下游查询决定是否采信。
                        log.warning("  入库但标记 flagged: %s (%s)",
                                    rec.iupac_short, rec.qc_notes)
                    gle.insert_record(conn, rec)
                    n += 1
                conn.close()
                log.info("  入库 %d 条（跳过 %d 条）", n, skipped)
            except Exception as e:  # noqa: BLE001
                log.error("  入库失败: %s", e)
                items[-1]["qc_issues"] = (items[-1]["qc_issues"] or "") + f" | DB_ERR: {e}"

    csv_path, md_path = dump_report(items, out_dir)
    log.info("完成: %d/%d 份有记录, 共 %d 条", 
             sum(1 for it in items if it["records"]), len(items),
             sum(it["records"] for it in items))
    print(f"\n报告: {md_path}\n     {csv_path}")


def asdict_min(res):
    return {"seq": res.residue_seq, "mono": res.monosaccharide_name,
            "ring": res.ring_form, "anomer": res.anomer,
            "reducing_end": res.is_reducing_end, "link": res.parent_carbon}


# ---------------------------------------------------------------------------
# 自检（用主脚本内置 SAMPLE_TEXT 模拟"PDF"跑批量全链路, 不连库）
# ---------------------------------------------------------------------------
def self_test():
    gle = load_etl(None)

    # monkeypatch 抽取函数: 每个假文件名对应独立样例文本
    orig = gle.extract_pdf_text
    SAMPLES = {
        "demo_oligo_glc.pdf": (
            "Compound 1: Glcα1-4Glc\n"
            "1H NMR (500 MHz, D2O) δ 5.41 (d, J = 3.8 Hz, 1H, H-1'), "
            "4.66 (d, J = 7.9 Hz, 1H, H-1β).\n"
            "13C NMR (126 MHz, D2O) δ 100.6 (C-1'), 96.8 (C-1), 78.3 (C-4).\n"
            "HSQC NMR (500 MHz, D2O) δ 5.41/100.6 (H-1', C-1').\n"
            "HMBC NMR (500 MHz, D2O) δ 5.41/78.3 (H-1', C-4).\n"),
        "demo_mono_glc.pdf": (
            "β-D-Glcp\n"
            "1H NMR (500 MHz, D2O) δ 4.64 (d, J = 7.9 Hz, 1H, H-1).\n"
            "13C NMR (126 MHz, D2O) δ 96.8 (C-1).\n"),
        "demo_poly_cellulose.pdf": (
            "→4)-β-D-Glcp-(1→4)-β-D-Glcp-(1→\n"
            "1H NMR (500 MHz, D2O) δ 4.52 (br s, 1H, H-1), 3.31 (m, 1H, H-2).\n"
            "13C NMR (126 MHz, D2O) δ 103.8 (C-1), 79.5 (C-4).\n"),
    }
    gle.extract_pdf_text = lambda p: SAMPLES[os.path.basename(p)]
    tmp = os.path.join(os.path.dirname(os.path.abspath(__file__)) or ".", "etl_reports")
    items = []
    fake = ["demo_oligo_glc.pdf", "demo_mono_glc.pdf", "demo_poly_cellulose.pdf"]
    for f in fake:
        _t, _b, recs = parse_one_pdf(gle, f, "10.xx", "J. Test", 2020)
        items.append({"file": f, "records": len(recs),
                      "flagged": sum(1 for r in recs if r.qc_status == "flagged"),
                      "passed": sum(1 for r in recs if r.qc_status == "passed"),
                      "types": ",".join(sorted({r.sugar_type for r in recs})),
                      "qc_issues": "; ".join("|".join(r.qc_notes) for r in recs if r.qc_notes),
                      "detail": ", ".join(r.iupac_short for r in recs)})
    csv_path, md_path = dump_report(items, tmp)
    assert items[0]["records"] == 1 and "oligo" in items[0]["types"], items[0]
    assert items[1]["records"] == 1 and "mono" in items[1]["types"], items[1]
    assert any("poly" in it["types"] for it in items), items
    gle.extract_pdf_text = orig
    print("SELF-TEST PASS")
    print(f"样例报告: {md_path}")
    print(f"         {csv_path}")


if __name__ == "__main__":
    main()
