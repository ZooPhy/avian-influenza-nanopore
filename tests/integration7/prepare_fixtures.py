#!/usr/bin/env python3
from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "tests" / "integration7" / "work"
PROJECT = WORK / "irma" / "project"


def write_fasta(path: Path, name: str, length: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    seq = ("ACGT" * ((length // 4) + 1))[:length]
    path.write_text(f">{name}\n{seq}\n", encoding="utf-8")


def write_fake_bam(path: Path, size: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"BAM_FIXTURE_" + (b"X" * size))


def write_coverage(path: Path, depths: list[int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["Position\tCoverage Depth"]
    lines.extend(f"{i+1}\t{depth}" for i, depth in enumerate(depths))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    if WORK.exists():
        shutil.rmtree(WORK)

    # HA has two complete candidates. A_HA_H5_HIGH should win because its
    # parseable coverage median is higher, despite the lower candidate being
    # lexicographically earlier.
    write_fasta(PROJECT / "assemblies" / "A_HA_H5_LOW.fasta", "A_HA_H5_LOW", 1650)
    write_fake_bam(PROJECT / "alignments" / "A_HA_H5_LOW.bam", 200)
    write_coverage(
        PROJECT / "coverage" / "A_HA_H5_LOW-coverage.txt",
        [20, 25, 30, 25, 20],
    )

    write_fasta(PROJECT / "assemblies" / "A_HA_H5_HIGH.fasta", "A_HA_H5_HIGH", 1700)
    write_fake_bam(PROJECT / "alignments" / "A_HA_H5_HIGH.bam", 250)
    write_coverage(
        PROJECT / "coverage" / "A_HA_H5_HIGH-coverage.txt",
        [90, 100, 110, 100, 95],
    )

    # NA is FASTA_ONLY.
    write_fasta(PROJECT / "assemblies" / "A_NA_N1.fasta", "A_NA_N1", 1400)

    # PB2 is BAM_ONLY.
    write_fake_bam(PROJECT / "alignments" / "A_PB2_TEST.bam", 150)

    # Remaining segments are absent.
    print(f"Prepared checkpoint-aware IRMA fixture under {PROJECT}")


if __name__ == "__main__":
    main()
