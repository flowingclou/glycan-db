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


# ---------------------------------------------------------------------------
# 报告生成
# ---------------------------------------------------------------------------
def dump_report(items, out_dir):
    """items: list of dict（每份 PDF 一行的汇总）。返回 (csv_path, md_path)。"""
    os.makedirs(out_dir, exist_ok=True)
    st = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_path = os.path.join(out_dir, f"manifest_{st}.csv")
    md_path = os.path.join(out_dir, f"manifest_{st}.md")
    fields = ["file", "records", "flagged", "passed",
              "types", "qc_issues", "detail"]
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for it in items:
            w.writerow(it)
    # Markdown 摘要
    ok = sum(1 for it in items if it["records"] > 0)
    total = sum(it["records"] for it in items)
    lines = [
        "# glycan_etl 批量执行报告",
        f"生成时间: {datetime.datetime.now().isoformat(timespec='seconds')}",
        "",
        f"- 处理 PDF 份数: **{len(items)}**",
        f"- 识别出记录的总 PDF 份数: **{ok}** / {len(items)}",
        f"- 共识别结构-谱图对: **{total}** 条",
        "",
        "| 文件 | 记录数 | flagged | passed | 类型分布 | QC问题 | 明细 |",
        "|------|-------:|-------:|-------:|---------|--------|------|",
    ]
    for it in items:
        lines.append(
            f"| {os.path.basename(it['file'])} | {it['records']} | {it['flagged']} | "
            f"{it['passed']} | {it['types'] or '-'} | {it['qc_issues'] or '-'} | "
            f"{(it['detail'] or '-')[:120]} |"
        )
    lines.append("")
    lines.append("> 0 条的 PDF 请优先检查文本层质量或结构写法；详细 JSON 见同目录 *_records.json")
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
    for i, pdf in enumerate(pdfs, 1):
        base = os.path.basename(pdf)
        m = meta.get(base) or {}
        doi = m.get("doi", default_doi)
        journal = m.get("journal", default_journal)
        year = m.get("year", default_year)
        log.info("[%d/%d] %s (doi=%s)", i, len(pdfs), base, doi)
        try:
            text, blocks, records = parse_one_pdf(gle, pdf, doi, journal, year)
        except Exception as e:  # noqa: BLE001
            log.error("解析失败: %s (%s)", base, e)
            items.append({"file": pdf, "records": 0, "flagged": 0, "passed": 0,
                          "types": "-", "qc_issues": f"ERROR: {e}", "detail": ""})
            continue
        flagged = sum(1 for r in records if r.qc_status == "flagged")
        passed = sum(1 for r in records if r.qc_status == "passed")
        types = ",".join(sorted({r.sugar_type for r in records}))
        issues = "; ".join(" | ".join(r.qc_notes) for r in records if r.qc_notes)
        detail = ", ".join(f"{r.iupac_short}({r.sugar_type},{len(r.experiments)}exp)"
                           for r in records)
        items.append({
            "file": pdf, "records": len(records), "flagged": flagged, "passed": passed,
            "types": types, "qc_issues": issues[:400], "detail": detail,
        })
        if args.save_json:
            with open(os.path.join(out_dir, f"{os.path.splitext(base)[0]}_records.json"),
                      "w", encoding="utf-8") as f:
                json.dump([{"iupac_short": r.iupac_short, "sugar_type": r.sugar_type,
                            "anomer": r.anomer, "n_experiments": len(r.experiments),
                            "n_2d_peaks": sum(len(e.peaks_2d) for e in r.experiments),
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
                n = 0
                for rec in records:
                    if rec.qc_status == "flagged":
                        continue
                    gle.insert_record(conn, rec)
                    n += 1
                conn.close()
                log.info("  入库 %d 条", n)
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
