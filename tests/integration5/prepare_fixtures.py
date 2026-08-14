#!/usr/bin/env python3
from __future__ import annotations

import gzip
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "tests" / "integration5" / "work"
MODEL = "dna_r10.4.1_e8.2_400bps_hac@v5.0.0"


def main() -> None:
    if WORK.exists():
        shutil.rmtree(WORK)

    (WORK / "data").mkdir(parents=True)
    (WORK / "coverage_flags").mkdir(parents=True)
    (WORK / "irma").mkdir(parents=True)

    with gzip.open(WORK / "data" / "smoke.fastq.gz", "wt", encoding="utf-8") as handle:
        for i in range(3):
            handle.write(
                f"@read{i+1} basecall_model_version_id={MODEL}\n"
                "ACGTACGTACGTACGT\n"
                "+\n"
                "FFFFFFFFFFFFFFFF\n"
            )

    (WORK / "coverage_flags" / "HA.flag").write_text("PASS\n", encoding="utf-8")
    (WORK / "coverage_flags" / "NA.flag").write_text("FAIL\n", encoding="utf-8")

    (WORK / "irma" / "HA.fasta").write_text(
        ">HA_IRMA\nACGTACGTACGTACGT\n",
        encoding="utf-8",
    )
    (WORK / "irma" / "NA.fasta").write_text(
        ">NA_IRMA\nTTTTCCCCAAAAGGGG\n",
        encoding="utf-8",
    )

    print(f"Prepared Phase 5 Medaka fixture under {WORK}")


if __name__ == "__main__":
    main()
