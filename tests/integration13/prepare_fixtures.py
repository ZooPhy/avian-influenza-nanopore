#!/usr/bin/env python3
from __future__ import annotations

import gzip
import random
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / "tests" / "integration13"
WORK = HERE / "work"
DATA = WORK / "data"
SOURCE = WORK / "source"

SAMPLE = "release_smoke"
MODEL = "dna_r10.4.1_e8.2_400bps_hac@v5.0.0"

SPECS = {
    "HA":  {"contig": "A_HA_H5_PHASE13", "length": 1700, "depth": 60},
    "NA":  {"contig": "A_NA_N1_PHASE13", "length": 1400, "depth": 20},
    "PB2": {"contig": "A_PB2_PHASE13",    "length": 2300, "depth": 20},
    "PB1": {"contig": "A_PB1_PHASE13",    "length": 2300, "depth": 20},
    "PA":  {"contig": "A_PA_PHASE13",     "length": 2200, "depth": 20},
    "NP":  {"contig": "A_NP_PHASE13",     "length": 1500, "depth": 20},
    "MP":  {"contig": "A_MP_PHASE13",     "length": 1000, "depth": 20},
    "NS":  {"contig": "A_NS_PHASE13",     "length":  850, "depth": 20},
}


def sequence_for(segment: str, length: int) -> str:
    rng = random.Random(1313 + sum(ord(c) for c in segment))
    return "".join(rng.choice("ACGT") for _ in range(length))


def write_sam(path: Path, contig: str, sequence: str, depth: int) -> None:
    quality = "I" * len(sequence)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("@HD\tVN:1.6\tSO:coordinate\n")
        handle.write(f"@SQ\tSN:{contig}\tLN:{len(sequence)}\n")
        for i in range(depth):
            handle.write(
                f"{segment_read_name(contig, i)}\t0\t{contig}\t1\t60\t"
                f"{len(sequence)}M\t*\t0\t0\t{sequence}\t{quality}\n"
            )


def segment_read_name(contig: str, index: int) -> str:
    return f"{contig.replace('_', '-')}-r{index + 1}"


def main() -> None:
    if WORK.exists():
        shutil.rmtree(WORK)

    DATA.mkdir(parents=True)
    SOURCE.mkdir(parents=True)

    sequences = {}
    for segment, spec in SPECS.items():
        sequence = sequence_for(segment, spec["length"])
        sequences[segment] = sequence

        (SOURCE / f"{segment}.fasta").write_text(
            f">{spec['contig']}\n{sequence}\n",
            encoding="utf-8",
        )
        write_sam(
            SOURCE / f"{segment}.sam",
            spec["contig"],
            sequence,
            spec["depth"],
        )

    # The raw Nanopore-style input is a tiny but valid FASTQ. It is deliberately
    # long enough and high quality to survive production Porechop/fastplong.
    ha = sequences["HA"]
    quality = "I" * len(ha)
    with gzip.open(DATA / f"{SAMPLE}.fastq.gz", "wt", encoding="utf-8") as handle:
        for i in range(60):
            handle.write(
                f"@raw{i+1} basecall_model_version_id={MODEL}\n"
                f"{ha}\n+\n{quality}\n"
            )

    (WORK / "metadata.tsv").write_text(
        "sample_id\thost_common_name\thost_species\tsample_type\t"
        "collection_date\tcollection_location\tstate\tcountry\tflyway\tproject_code\n"
        "release_smoke\tMallard\tAnas platyrhynchos\toropharyngeal swab\t"
        "2026-08-14\tPhoenix, Arizona\tArizona\tUSA\tPacific\tP13-E2E\n",
        encoding="utf-8",
    )

    # HA is the only sequence in the tiny BLAST reference database.
    (SOURCE / "blast_reference.fasta").write_text(
        ">PHASE13_H5_REFERENCE synthetic_H5_phase13_reference\n"
        + ha
        + "\n",
        encoding="utf-8",
    )

    print(
        "Prepared Phase 13 raw FASTQ, metadata, segment fixtures, and BLAST reference."
    )


if __name__ == "__main__":
    main()
