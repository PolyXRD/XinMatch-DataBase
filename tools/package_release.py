"""Package built ``.xinDB`` libraries into release zips + SHA256SUMS.txt.

Default: every ``data/*.xinDB``.  Deflate level 1 keeps ~1.2x on the mostly
float32 payload at ~35 MB/s, so a five-library set takes a few minutes.

Usage::

    python -m tools.package_release --out _release/20261104
    python -m tools.package_release --out _release/20261104 COD-Inorg-20261020.xinDB
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import time
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="输出目录，如 _release/20261104")
    ap.add_argument("--data-dir", default=str(DATA), help="库目录（默认 <repo>/data）")
    ap.add_argument("names", nargs="*", help="要打包的文件名（默认 data 下全部 *.xinDB）")
    args = ap.parse_args(argv[1:])

    data_dir = Path(args.data_dir)
    out_dir = Path(args.out)
    if not out_dir.is_absolute():
        out_dir = (REPO / out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    names = args.names or sorted(p.name for p in data_dir.glob("*.xinDB"))
    # ICDD PDF-2 是商品库，**永不进 Release**。默认 glob 会把不小心放进 data/
    # 的 PDF2 库一起打进去，所以在默认路径上直接剔除；显式点名则说明使用者
    # 知道自己在做什么，不再拦。
    if not args.names:
        licensed = [n for n in names if "PDF2" in n.upper()]
        if licensed:
            print(
                f"[SKIP] 跳过 ICDD 授权库（不得发布）：{', '.join(licensed)}",
                file=sys.stderr,
            )
            names = [n for n in names if n not in licensed]
    if not names:
        print(f"[ERR] no .xinDB under {data_dir}", file=sys.stderr)
        return 1

    sums: list[tuple[str, str]] = []
    for name in names:
        src = data_dir / name
        if not src.exists():
            print(f"[MISS] {src}", file=sys.stderr)
            return 1
        dst = out_dir / f"{name}.zip"
        print(f"[ZIP ] {name} ({src.stat().st_size:,} B) -> {dst}", flush=True)
        t0 = time.perf_counter()
        with zipfile.ZipFile(dst, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=1) as zf:
            zf.write(src, arcname=name)
        digest = sha256(dst)
        sums.append((digest, dst.name))
        print(
            f"[OK  ] {dst.name}: {dst.stat().st_size:,} B  ratio={src.stat().st_size / dst.stat().st_size:.2f}x  "
            f"{src.stat().st_size / (time.perf_counter() - t0) / 1e6:.0f} MB/s  sha256={digest}",
            flush=True,
        )

    sums_path = out_dir / "SHA256SUMS.txt"
    sums_path.write_text("\n".join(f"{h}  {n}" for h, n in sums) + "\n", encoding="utf-8")
    print(f"[SUMS] {sums_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
