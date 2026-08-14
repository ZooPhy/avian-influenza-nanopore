#!/usr/bin/env python3
from __future__ import annotations

import gzip
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "tests" / "integration8" / "work"
SOURCE = WORK / "source"
DATA = WORK / "data"
SAMPLE = "qc_checkpoint"

SPECS = {
    "HA":  {"contig": "A_HA_H5",   "length": 1700, "depth": 60, "n_count": 0},
    "NA":  {"contig": "A_NA_N1",   "length": 1400, "depth": 20, "n_count": 0},
    "PB2": {"contig": "A_PB2_TEST", "length": 2450, "depth": 60, "n_count": 0},
    "PB1": {"contig": "A_PB1_TEST", "length": 2100, "depth": 60, "n_count": 0},
    "PA":  {"contig": "A_PA_TEST",  "length": 2200, "depth": 60, "n_count": 50},
    "NP":  {"contig": "A_NP_TEST",  "length": 1500, "depth": 60, "n_count": 0},
    "MP":  {"contig": "A_MP_TEST",  "length": 1000, "depth": 60, "n_count": 0},
    "NS":  {"contig": "A_NS_TEST",  "length":  850, "depth": 60, "n_count": 0},
}


def consensus(length: int, n_count: int) -> str:
    if n_count > length:
        raise ValueError("n_count exceeds sequence length")
    return ("N" * n_count) + ("A" * (length - n_count))


def write_sam(path: Path, contig: str, length: int, depth: int) -> None:
    read_seq = "A" * length
    quality = "I" * length
    with path.open("w", encoding="utf-8") as handle:
        handle.write("@HD\tVN:1.6\tSO:coordinate\n")
        handle.write(f"@SQ\tSN:{contig}\tLN:{length}\n")
        for i in range(depth):
            handle.write(
                f"r{i+1}\t0\t{contig}\t1\t60\t{length}M\t*\t0\t0\t"
                f"{read_seq}\t{quality}\n"
            )


def main() -> None:
    if WORK.exists():
        shutil.rmtree(WORK)
    SOURCE.mkdir(parents=True)
    DATA.mkdir(parents=True)

    for segment, spec in SPECS.items():
        contig = spec["contig"]
        length = spec["length"]
        seq = consensus(length, spec["n_count"])

        (SOURCE / f"{segment}.fasta").write_text(
            f">{contig}\n{seq}\n",
            encoding="utf-8",
        )
        write_sam(
            SOURCE / f"{segment}.sam",
            contig,
            length,
            spec["depth"],
        )

    with gzip.open(DATA / f"{SAMPLE}.fastq.gz", "wt", encoding="utf-8") as handle:
        handle.write(
            "@fixture basecall_model_version_id="
            "dna_r10.4.1_e8.2_400bps_hac@v5.0.0\n"
            "ACGT\n+\nIIII\n"
        )

    (WORK / "metadata.tsv").write_text(
        "sample_id\nqc_checkpoint\n",
        encoding="utf-8",
    )

    print(f"Prepared Phase 8 source fixtures under {WORK}")


if __name__ == "__main__":
    main()
