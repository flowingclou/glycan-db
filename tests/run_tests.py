# -*- coding: utf-8 -*-
"""glycan-db 仓库自检入口。

功能:
  1) 递归 py_compile 校验仓库内所有 .py 脚本;
  2) 运行 `pipelines/batch_etl.py --self-test` 验证解析引擎与批量管道可用。

用法:
  python3 tests/run_tests.py
"""
import pathlib
import py_compile
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent


def compile_all():
    py_files = sorted(ROOT.rglob("*.py"))
    if not py_files:
        raise SystemExit("未找到任何 .py 文件")
    for f in py_files:
        py_compile.compile(str(f), doraise=True)
    return py_files


def run_self_test():
    script = ROOT / "pipelines" / "batch_etl.py"
    proc = subprocess.run(
        [sys.executable, str(script), "--self-test"],
        capture_output=True,
        text=True,
    )
    print(proc.stdout)
    if proc.returncode != 0:
        print(proc.stderr, file=sys.stderr)
        return False
    return True


def main() -> int:
    print("== 1/2  py_compile 全仓校验 ==")
    py_files = compile_all()
    print(f"OK: {len(py_files)} 个 .py 文件全部编译通过\n")

    print("== 2/2  batch_etl --self-test ==")
    if not run_self_test():
        print("自检失败", file=sys.stderr)
        return 1
    print("自检通过")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001
        print(f"[fatal] {type(exc).__name__}: {exc}", file=sys.stderr)
        sys.exit(1)
