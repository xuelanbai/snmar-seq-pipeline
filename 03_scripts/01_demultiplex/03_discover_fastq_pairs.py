#!/usr/bin/env python3
"""Step 03: discover paired FASTQ files below the four multipath raw-data folders."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


SUFFIX_PAIRS = (
    ("_1.fq.gz", "_2.fq.gz"),
    ("_R1.fq.gz", "_R2.fq.gz"),
    ("_R1.fastq.gz", "_R2.fastq.gz"),
    ("_1.fastq.gz", "_2.fastq.gz"),
)


def paired_path(r1: Path) -> Path | None:
    text = str(r1)
    for left, right in SUFFIX_PAIRS:
        if text.endswith(left):
            return Path(text[: -len(left)] + right)
    return None


def is_r1(path: Path) -> bool:
    return any(str(path).endswith(left) for left, _ in SUFFIX_PAIRS)


def safe_pair_id(relative_r1: Path) -> str:
    stem = relative_r1.name
    for left, _ in SUFFIX_PAIRS:
        if stem.endswith(left):
            stem = stem[: -len(left)]
            break
    digest = hashlib.sha1(str(relative_r1).encode("utf-8")).hexdigest()[:10]
    cleaned = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in stem)
    return f"{cleaned}.{digest}"


def sample_name(path: Path, raw_root: Path) -> str:
    rel = path.relative_to(raw_root)
    parts = rel.parts
    if len(parts) >= 2:
        return parts[-2]
    return path.parent.name


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--undetermined-out", required=True, type=Path)
    args = parser.parse_args()

    raw_root = args.raw_root.resolve()
    if not raw_root.is_dir():
        raise SystemExit(f"[ERROR] raw root does not exist: {raw_root}")

    paired_rows: list[tuple[str, str, str, str]] = []
    undetermined_rows: list[tuple[str, str, str, str]] = []
    missing: list[str] = []

    for r1 in sorted(p for p in raw_root.rglob("*.gz") if is_r1(p)):
        r2 = paired_path(r1)
        assert r2 is not None
        if not r2.is_file():
            missing.append(str(r2))
            continue
        rel = r1.relative_to(raw_root)
        row = (safe_pair_id(rel), sample_name(r1, raw_root), str(r1), str(r2))
        if any(part.lower() == "undetermined" for part in rel.parts):
            undetermined_rows.append(row)
        else:
            paired_rows.append(row)

    if missing:
        preview = "\n".join(missing[:20])
        raise SystemExit(f"[ERROR] missing {len(missing)} R2 files; first entries:\n{preview}")
    if not paired_rows:
        raise SystemExit(f"[ERROR] no paired FASTQ files found below {raw_root}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    header = "pair_id\tsample\tR1\tR2\n"
    args.out.write_text(
        header + "".join("\t".join(row) + "\n" for row in paired_rows),
        encoding="utf-8",
    )
    args.undetermined_out.write_text(
        header + "".join("\t".join(row) + "\n" for row in undetermined_rows),
        encoding="utf-8",
    )
    print(f"[INFO] regular FASTQ pairs: {len(paired_rows)}")
    print(f"[INFO] Undetermined FASTQ pairs: {len(undetermined_rows)}")
    print(f"[INFO] pair table: {args.out}")


if __name__ == "__main__":
    main()
