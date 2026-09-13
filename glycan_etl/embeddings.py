#!/usr/bin/env python3
"""
glycan_embed_bge3.py — 使用硅基流动 SiliconFlow BAAI/bge-m3 为 nmr_shifts_1d 生成 1024 维 embedding 并写入库。
幂等可续跑：已有 embedding 的行跳过。失败自动重试(指数退避)。
用法:
  SILICONFLOW_API_KEY=sk-xxx python3 glycan_embed_bge3.py [--batch 20] [--limit 0]
"""
import os, sys, time, json, argparse, urllib.request, urllib.error

try:
    import psycopg2
except ImportError:
    sys.exit("缺少依赖 psycopg2，请先 pip install psycopg2-binary")

API_URL = "https://api.siliconflow.cn/v1/embeddings"
MODEL = "BAAI/bge-m3"
KEY = os.environ.get("SILICONFLOW_API_KEY", "")
DB = dict(host="localhost", port=5432, dbname="glycan_db", user=os.getenv("PGUSER", os.getenv("USER")))


def build_text(row):
    """为一条位移记录构造可检索语义文本(融合糖链上下文)。"""
    parts = [
        f"{row['nucleus']} NMR chemical shift {float(row['shift_ppm']):.3f} ppm",
    ]
    if row["assignment_position"]:
        parts.append(f"assigned to position {row['assignment_position']}")
    if row["multiplicity"]:
        parts.append(f"multiplicity {row['multiplicity']}")
    if row["j_coupling_hz"]:
        parts.append(f"J coupling {float(row['j_coupling_hz']):.2f} Hz")
    if row["is_anomeric"]:
        parts.append("anomeric proton/carbon")
    if row["iupac_short"]:
        parts.append(f"context: {row['iupac_short']}")
    if row["solvent"]:
        parts.append(f"solvent {row['solvent']}")
    return ", ".join(parts)


def call_embedding(texts):
    payload = {"model": MODEL, "input": texts}
    req = urllib.request.Request(
        API_URL,
        data=json.dumps(payload).encode(),
        headers={"Authorization": f"Bearer {KEY}", "Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=60) as r:
        body = json.loads(r.read())
    return [d["embedding"] for d in sorted(body["data"], key=lambda x: x["index"])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=int, default=20)
    ap.add_argument("--limit", type=int, default=0, help="0=全部")
    args = ap.parse_args()
    if not KEY:
        sys.exit("未设置环境变量 SILICONFLOW_API_KEY")

    conn = psycopg2.connect(**DB)
    conn.autocommit = False
    cur = conn.cursor()

    # 未向量化的目标行，join 糖链上下文
    sql = """
        SELECT s.shift_id, s.nucleus, s.shift_ppm, s.multiplicity,
               s.j_coupling_hz, s.assignment_position, s.is_anomeric,
               g.iupac_short, e.solvent
        FROM nmr_shifts_1d s
        JOIN nmr_experiments e ON e.experiment_id = s.experiment_id
        LEFT JOIN sugars g ON g.sugar_id = e.sugar_id
        WHERE s.embedding IS NULL
        ORDER BY s.shift_id
    """
    if args.limit:
        sql += f" LIMIT {args.limit}"
    cur.execute(sql)
    rows = cur.fetchall()
    cols = [d[0] for d in cur.description]
    total = len(rows)
    print(f"[1/2] 待向量化 {total} 条", flush=True)

    done, fail = 0, 0
    for i in range(0, total, args.batch):
        chunk = rows[i : i + args.batch]
        texts = [build_text(dict(zip(cols, r))) for r in chunk]
        emb = None
        for attempt in range(4):
            try:
                emb = call_embedding(texts)
                break
            except Exception as ex:
                wait = 2 ** attempt * 2
                print(f"  批 {i//args.batch} 第{attempt+1}次失败: {ex}，{wait}s 后重试", flush=True)
                time.sleep(wait)
        if emb is None:
            print(f"  批 {i//args.batch} 重试耗尽，跳过", flush=True)
            fail += len(chunk)
            continue
        cur.executemany(
            "UPDATE nmr_shifts_1d SET embedding=%s::vector WHERE shift_id=%s",
            [(e, r[0]) for e, r in zip(emb, chunk)],
        )
        conn.commit()
        done += len(emb)
        print(f"  已写入 {done}/{total} (失败 {fail})", flush=True)
        if (i // args.batch) % 5 == 4:
            time.sleep(1)  # 轻柔限流

    cur.close()
    conn.close()
    print(f"[2/2] 完成: 写入 {done}，失败 {fail}", flush=True)


if __name__ == "__main__":
    main()
