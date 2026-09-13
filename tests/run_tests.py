# -*- coding: utf-8 -*-
"""glycan-db 仓库自检入口。

功能:
  1) 递归 py_compile 校验仓库内所有 .py 脚本;
  2) 运行 `tests/test_regressions.py` 回归测试（防回退，不依赖数据库）;
  3) 运行 `pipelines/batch_etl.py --self-test` 验证解析引擎与批量管道可用。

用法:
  python3 tests/run_tests.py
"""
import pathlib
import py_compile
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent

# 仓库自身源码之外的内容不参与编译校验：虚拟环境里会混入成百上千个第三方
# .py 文件，既拖慢自检，也可能因第三方语法/依赖导致误报。
SKIP_DIRS = {".git", ".venv", "venv", "env", "ENV", "__pycache__",
             "build", "dist", ".eggs", "node_modules", "site-packages"}


def compile_all():
    py_files = sorted(
        p for p in ROOT.rglob("*.py")
        if not SKIP_DIRS.intersection(p.relative_to(ROOT).parts)
    )
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


def run_regressions():
    script = ROOT / "tests" / "test_regressions.py"
    proc = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True,
        text=True,
    )
    print(proc.stdout)
    if proc.returncode != 0:
        print(proc.stderr, file=sys.stderr)
        return False
    return True


def main() -> int:
    print("== 1/3  py_compile 全仓校验 ==")
    py_files = compile_all()
    print(f"OK: {len(py_files)} 个 .py 文件全部编译通过\n")

    print("== 2/3  回归测试（防回退） ==")
    if not run_regressions():
        print("回归测试失败", file=sys.stderr)
        return 1
    print("回归测试通过\n")

    print("== 3/3  batch_etl --self-test ==")
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
