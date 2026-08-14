#!/usr/bin/env python3
from __future__ import annotations

import gzip
import random
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "tests" / "integration9" / "work"
SOURCE = WORK / "source"
DATA = WORK / "data"

SAMPLE = "phase9"
CONTIG = "A_HA_H5_PHASE9"
MODEL = "dna_r10.4.1_e8.2_400bps_hac@v5.0.0"
LENGTH = 1700
DEPTH = 60


def make_reference() -> str:
    rng = random.Random(909)
    return "".join(rng.choice("ACGT") for _ in range(LENGTH))


def main() -> None:
    if WORK.exists():
        shutil.rmtree(WORK)

    SOURCE.mkdir(parents=True)
    DATA.mkdir(parents=True)

    reference = make_reference()

    (SOURCE / "HA.fasta").write_text(
        f">{CONTIG}\n{reference}\n",
        encoding="utf-8",
    )

    # Exact reference is also the only sequence in the tiny BLAST database.
    (SOURCE / "blast_reference.fasta").write_text(
        ">PHASE9_H5_REFERENCE synthetic_H5_phase9_reference\n"
        + reference
        + "\n",
        encoding="utf-8",
    )

    quality = "I" * LENGTH
    with (SOURCE / "HA.sam").open("w", encoding="utf-8") as handle:
        handle.write("@HD\tVN:1.6\tSO:coordinate\n")
        handle.write(f"@SQ\tSN:{CONTIG}\tLN:{LENGTH}\n")
        for i in range(DEPTH):
            handle.write(
                f"r{i+1}\t0\t{CONTIG}\t1\t60\t{LENGTH}M\t*\t0\t0\t"
                f"{reference}\t{quality}\n"
            )

    # Production model resolution reads the original sample FASTQ. The reads are
    # identical to the BAM fixture so Medaka receives a coherent tiny dataset.
    with gzip.open(DATA / f"{SAMPLE}.fastq.gz", "wt", encoding="utf-8") as handle:
        for i in range(DEPTH):
            handle.write(
                f"@r{i+1} basecall_model_version_id={MODEL}\n"
                f"{reference}\n"
                "+\n"
                f"{quality}\n"
            )

    (WORK / "metadata.tsv").write_text(
        "sample_id\nphase9\n",
        encoding="utf-8",
    )

    print(
        f"Prepared Phase 9 fixture: one {LENGTH}-nt HA segment at "
        f"{DEPTH}x nominal depth under {WORK}"
    )


if __name__ == "__main__":
    main()
