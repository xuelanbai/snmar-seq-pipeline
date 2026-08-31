#!/usr/bin/env python3
"""Step 03: discover RNA R1/R2 FASTQ pairs and validate samples 1..96."""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path


R1_RE = re.compile(
    r"^(?P<prefix>.+)_R1(?P<tail>(?:_[^.]+)?)\.(?P<ext>f(?:ast)?q\.gz)$",
    re.IGNORECASE,
)
SAMPLE_LANE_RE = re.compile(r"^(?P<sample>.+?)(?:_S\d+)?_L\d+$", re.IGNORECASE)
SAMPLE_NUMBER_RE = re.compile(r"-(\d+)$")


def paired_path(r1: Path) -> Path | None:
    match = R1_RE.match(r1.name)
    if match is None:
        return None
    r2_name = (
        f"{match.group('prefix')}_R2{match.group('tail')}."
        f"{match.group('ext')}"
    )
    return r1.with_name(r2_name)


def sample_name(r1: Path) -> str:
    match = R1_RE.match(r1.name)
    assert match is not None
    prefix = match.group("prefix")
    lane_match = SAMPLE_LANE_RE.match(prefix)
    return lane_match.group("sample") if lane_match else prefix


def safe_pair_id(relative_r1: Path) -> str:
    match = R1_RE.match(relative_r1.name)
    assert match is not None
    stem = f"{match.group('prefix')}{match.group('tail')}"
    cleaned = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in stem)
    digest = hashlib.sha1(str(relative_r1).encode("utf-8")).hexdigest()[:10]
    return f"{cleaned}.{digest}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--expected-samples", type=int, default=96)
    args = parser.parse_args()

    raw_root = args.raw_root.resolve()
    if not raw_root.is_dir():
        raise SystemExit(f"[ERROR] RNA raw root does not exist: {raw_root}")

    rows: list[tuple[str, str, str, str]] = []
    missing: list[str] = []
    for r1 in sorted(raw_root.rglob("*.gz")):
        r2 = paired_path(r1)
        if r2 is None:
            continue
        if not r2.is_file():
            missing.append(str(r2))
            continue
        relative_r1 = r1.relative_to(raw_root)
        rows.append(
            (
                safe_pair_id(relative_r1),
                sample_name(r1),
                str(r1),
                str(r2),
            )
        )

    if missing:
        preview = "\n".join(missing[:20])
        raise SystemExit(
            f"[ERROR] missing {len(missing)} RNA R2 file(s); first entries:\n{preview}"
        )
    if not rows:
        raise SystemExit(f"[ERROR] no RNA FASTQ pairs found below {raw_root}")

    samples = sorted({row[1] for row in rows})
    observed_numbers: set[int] = set()
    unnumbered: list[str] = []
    for sample in samples:
        match = SAMPLE_NUMBER_RE.search(sample)
        if match is None:
            unnumbered.append(sample)
        else:
            observed_numbers.add(int(match.group(1)))

    expected_numbers = set(range(1, args.expected_samples + 1))
    missing_numbers = sorted(expected_numbers - observed_numbers)
    unexpected_numbers = sorted(observed_numbers - expected_numbers)
    if unnumbered or missing_numbers or unexpected_numbers:
        raise SystemExit(
            "[ERROR] RNA sample-number validation failed: "
            f"unnumbered={unnumbered[:10]}, "
            f"missing={missing_numbers}, unexpected={unexpected_numbers}"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    header = "pair_id\tsample\tR1\tR2\n"
    args.out.write_text(
        header + "".join("\t".join(row) + "\n" for row in rows),
        encoding="utf-8",
    )
    print(f"[INFO] RNA FASTQ pairs: {len(rows)}")
    print(f"[INFO] RNA samples: {len(samples)} (validated 1..{args.expected_samples})")
    print(f"[INFO] pair table: {args.out}")


if __name__ == "__main__":
    main()

