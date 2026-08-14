#!/usr/bin/env python3
from __future__ import annotations

import random
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "tests" / "integration4" / "work"


def dna(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def main() -> None:
    if WORK.exists():
        shutil.rmtree(WORK)

    (WORK / "queries").mkdir(parents=True)
    (WORK / "flags").mkdir(parents=True)

    ha = dna(1700, 42)
    decoy = dna(1700, 99)
    na = dna(1400, 123)

    (WORK / "reference.fasta").write_text(
        ">TEST_HA_REF H5 synthetic exact-match reference\n"
        + ha
        + "\n>DECOY_REF synthetic decoy reference\n"
        + decoy
        + "\n",
        encoding="utf-8",
    )

    (WORK / "queries" / "HA.fasta").write_text(
        ">QUERY_HA\n" + ha + "\n",
        encoding="utf-8",
    )
    (WORK / "queries" / "NA.fasta").write_text(
        ">QUERY_NA\n" + na + "\n",
        encoding="utf-8",
    )

    (WORK / "flags" / "HA.flag").write_text("PASS\n", encoding="utf-8")
    (WORK / "flags" / "NA.flag").write_text("FAIL\n", encoding="utf-8")

    print(f"Prepared Phase 4 BLAST fixture under {WORK}")


if __name__ == "__main__":
    main()
