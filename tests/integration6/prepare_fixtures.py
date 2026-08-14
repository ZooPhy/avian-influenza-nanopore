#!/usr/bin/env python3
from __future__ import annotations

import gzip
import random
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "tests" / "integration6" / "work"
MODEL = "dna_r10.4.1_e8.2_400bps_hac@v5.0.0"


def make_reference(length: int = 1400) -> str:
    rng = random.Random(606)
    return "".join(rng.choice("ACGT") for _ in range(length))


def mutate(seq: str, rng: random.Random, rate: float = 0.012) -> str:
    bases = "ACGT"
    out = []
    for base in seq:
        if rng.random() < rate:
            choices = [b for b in bases if b != base]
            out.append(rng.choice(choices))
        else:
            out.append(base)
    return "".join(out)


def main() -> None:
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True)

    reference = make_reference()
    (WORK / "HA.fasta").write_text(
        ">HA_SYNTHETIC\n" + reference + "\n",
        encoding="utf-8",
    )

    rng = random.Random(607)
    with gzip.open(WORK / "reads.fastq.gz", "wt", encoding="utf-8") as handle:
        for i in range(40):
            read = mutate(reference, rng)
            handle.write(
                f"@read{i+1} basecall_model_version_id={MODEL}\n"
                f"{read}\n"
                "+\n"
                f"{'I' * len(read)}\n"
            )

    print(
        f"Prepared Phase 6 fixture: 40 x {len(reference)}-nt reads "
        f"against one synthetic HA contig under {WORK}"
    )


if __name__ == "__main__":
    main()
