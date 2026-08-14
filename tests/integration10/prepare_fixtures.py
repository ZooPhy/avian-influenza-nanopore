#!/usr/bin/env python3
from __future__ import annotations

import gzip
import json
import random
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "tests" / "integration10" / "work"
SOURCE = WORK / "source"
DATA = WORK / "data"
RESULTS = WORK / "results"
SAMPLE = "phase10"
MODEL = "dna_r10.4.1_e8.2_400bps_hac@v5.0.0"

SPECS = {
    "HA":  {"contig": "A_HA_H5_PHASE10", "length": 1700, "depth": 60},
    "NA":  {"contig": "A_NA_N1_PHASE10", "length": 1400, "depth": 20},
    "PB2": {"contig": "A_PB2_PHASE10",    "length": 2300, "depth": 20},
    "PB1": {"contig": "A_PB1_PHASE10",    "length": 2300, "depth": 20},
    "PA":  {"contig": "A_PA_PHASE10",     "length": 2200, "depth": 20},
    "NP":  {"contig": "A_NP_PHASE10",     "length": 1500, "depth": 20},
    "MP":  {"contig": "A_MP_PHASE10",     "length": 1000, "depth": 20},
    "NS":  {"contig": "A_NS_PHASE10",     "length":  850, "depth": 20},
}


def deterministic_sequence(segment: str, length: int) -> str:
    rng = random.Random(1010 + sum(ord(c) for c in segment))
    return "".join(rng.choice("ACGT") for _ in range(length))


def write_sam(path: Path, contig: str, sequence: str, depth: int) -> None:
    quality = "I" * len(sequence)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("@HD\tVN:1.6\tSO:coordinate\n")
        handle.write(f"@SQ\tSN:{contig}\tLN:{len(sequence)}\n")
        for i in range(depth):
            handle.write(
                f"r{i+1}\t0\t{contig}\t1\t60\t{len(sequence)}M\t*\t0\t0\t"
                f"{sequence}\t{quality}\n"
            )


def main() -> None:
    if WORK.exists():
        shutil.rmtree(WORK)

    SOURCE.mkdir(parents=True)
    DATA.mkdir(parents=True)

    sequences = {}
    for segment, spec in SPECS.items():
        sequence = deterministic_sequence(segment, spec["length"])
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

    # Only HA can pass the 50x QC threshold. It is also the exact BLAST reference.
    (SOURCE / "blast_reference.fasta").write_text(
        ">PHASE10_H5_REFERENCE synthetic_H5_phase10_reference\n"
        + sequences["HA"]
        + "\n",
        encoding="utf-8",
    )

    # Production model resolution reads the original sample FASTQ.
    with gzip.open(DATA / f"{SAMPLE}.fastq.gz", "wt", encoding="utf-8") as handle:
        quality = "I" * len(sequences["HA"])
        for i in range(60):
            handle.write(
                f"@ha_read{i+1} basecall_model_version_id={MODEL}\n"
                f"{sequences['HA']}\n+\n{quality}\n"
            )

    # Top-level metadata file satisfies production sample discovery/configuration.
    (WORK / "metadata.tsv").write_text(
        "sample_id\thost_species\tcollection_location\tflyway\tproject_code\tsegments_pass\n"
        "phase10\tAnas platyrhynchos\tPhoenix, Arizona\tPacific\tP10-REPORT\tfield_value\n",
        encoding="utf-8",
    )

    sample_dir = RESULTS / SAMPLE

    # Direct sample_summary prerequisites that are upstream of the reporting
    # boundary but are not the focus of Phase 10.
    metadata_dir = sample_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    (metadata_dir / f"{SAMPLE}.metadata.tsv").write_text(
        "sample_id\thost_species\tcollection_location\tflyway\tproject_code\tsegments_pass\n"
        "phase10\tAnas platyrhynchos\tPhoenix, Arizona\tPacific\tP10-REPORT\tfield_value\n",
        encoding="utf-8",
    )

    fastplong_dir = sample_dir / "fastplong"
    fastplong_dir.mkdir(parents=True, exist_ok=True)
    fastplong = {
        "summary": {
            "before_filtering": {
                "total_reads": 1000,
                "total_bases": 1500000,
                "q20_rate": 0.91,
                "q30_rate": 0.72,
            },
            "after_filtering": {
                "total_reads": 900,
                "total_bases": 1400000,
                "q20_rate": 0.95,
                "q30_rate": 0.79,
            },
        },
        "filtering_result": {
            "passed_filter_reads": 900,
            "low_quality_reads": 60,
            "too_short_reads": 40,
        },
        "adapter_cutting": {
            "adapter_trimmed_reads": 0,
        },
    }
    (fastplong_dir / "report.json").write_text(
        json.dumps(fastplong),
        encoding="utf-8",
    )

    genoflu_dir = sample_dir / "genoflu"
    genoflu_dir.mkdir(parents=True, exist_ok=True)
    (genoflu_dir / "GenoFLU.tsv").write_text(
        "sample\tstatus\nphase10\tSKIPPED_FIXTURE\n",
        encoding="utf-8",
    )

    vadr_dir = sample_dir / "vadr"
    vadr_dir.mkdir(parents=True, exist_ok=True)
    (vadr_dir / f"{SAMPLE}.vadr.log").write_text(
        "Phase 10 reporting fixture: VADR not executed.\n",
        encoding="utf-8",
    )

    # concat_consensus requires a VCF for every normalized READY segment.
    # Phase 10 is not a VCF test, so explicit fixture files prevent production
    # medaka_vcf from being scheduled while still allowing real consensus
    # generation to flow into concat_consensus.
    for segment in SPECS:
        medaka_dir = sample_dir / "medaka" / segment
        medaka_dir.mkdir(parents=True, exist_ok=True)
        (medaka_dir / "variants.vcf").write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            encoding="utf-8",
        )
        (medaka_dir / "variants.status.tsv").write_text(
            "status\treason\nSKIPPED_FIXTURE\tphase10_reporting_boundary\n",
            encoding="utf-8",
        )

    print(
        "Prepared Phase 10 fixtures: eight READY IRMA-like segments, "
        "only HA QC-passing, plus reporting-boundary prerequisites."
    )


if __name__ == "__main__":
    main()
