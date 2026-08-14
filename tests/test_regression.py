from __future__ import annotations

import csv
import importlib.util
import json
import runpy
from pathlib import Path
from types import SimpleNamespace


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relative_path: str):
    path = REPO_ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def read_one_tsv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def test_check_coverage_warning_and_multiple_candidates_still_pass(tmp_path, monkeypatch):
    mod = load_module("check_coverage", "scripts/check_coverage.py")

    segments_dir = tmp_path / "segments"
    ha_dir = segments_dir / "HA"
    ha_dir.mkdir(parents=True)
    (ha_dir / "alignment.bam").write_bytes(b"placeholder")
    consensus = ha_dir / "consensus.fasta"
    consensus.write_text(">A_HA_H5\n" + "A" * 1909 + "\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "segment\tstatus\tcontig\tcandidate_count\tselection_status\n"
        "HA\tREADY\tA_HA_H5\t2\tMULTIPLE_CANDIDATES\n",
        encoding="utf-8",
    )

    flag = tmp_path / "HA.flag"
    stats = tmp_path / "HA.tsv"

    monkeypatch.setattr(mod, "samtools_depths", lambda _: ([100] * 1909, "A_HA_H5"))
    mod.snakemake = SimpleNamespace(
        input=SimpleNamespace(
            segments_dir=str(segments_dir), manifest=str(manifest), consensus=str(consensus)
        ),
        output=SimpleNamespace(flag=str(flag), stats=str(stats)),
        wildcards=SimpleNamespace(segment="HA"),
        params=SimpleNamespace(
            min_median_depth=50,
            min_breadth=0.95,
            expected_length_min=1600,
            expected_length_max=1800,
            max_n_fraction=0.01,
        ),
    )

    mod.main()
    row = read_one_tsv(stats)

    assert flag.read_text().strip() == "PASS"
    assert row["coverage_status"] == "PASS"
    assert row["length_status"] == "WARNING"
    assert row["n_content_status"] == "PASS"
    assert row["status"] == "PASS"
    assert row["n_candidates"] == "2"
    assert row["selection_status"] == "MULTIPLE_CANDIDATES"


def test_check_coverage_low_breadth_fails(tmp_path, monkeypatch):
    mod = load_module("check_coverage_low_breadth", "scripts/check_coverage.py")

    segments_dir = tmp_path / "segments"
    na_dir = segments_dir / "NA"
    na_dir.mkdir(parents=True)
    (na_dir / "alignment.bam").write_bytes(b"placeholder")
    consensus = na_dir / "consensus.fasta"
    consensus.write_text(">A_NA_N1\n" + "A" * 1350 + "\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "segment\tstatus\tcontig\tcandidate_count\tselection_status\n"
        "NA\tREADY\tA_NA_N1\t1\tUNIQUE\n",
        encoding="utf-8",
    )

    flag = tmp_path / "NA.flag"
    stats = tmp_path / "NA.tsv"
    depths = [100] * 1200 + [0] * 150  # 88.9% breadth at 50x
    monkeypatch.setattr(mod, "samtools_depths", lambda _: (depths, "A_NA_N1"))
    mod.snakemake = SimpleNamespace(
        input=SimpleNamespace(
            segments_dir=str(segments_dir), manifest=str(manifest), consensus=str(consensus)
        ),
        output=SimpleNamespace(flag=str(flag), stats=str(stats)),
        wildcards=SimpleNamespace(segment="NA"),
        params=SimpleNamespace(
            min_median_depth=50,
            min_breadth=0.95,
            expected_length_min=1300,
            expected_length_max=1500,
            max_n_fraction=0.01,
        ),
    )

    mod.main()
    row = read_one_tsv(stats)

    assert flag.read_text().strip() == "FAIL"
    assert row["coverage_status"] == "FAIL"
    assert row["length_status"] == "PASS"
    assert row["status"] == "FAIL"


def test_blast_summary_four_states(tmp_path):
    blast_dir = tmp_path / "blast"
    flags_dir = tmp_path / "flags"
    blast_dir.mkdir()
    flags_dir.mkdir()

    # HA: high confidence
    (flags_dir / "HA.flag").write_text("PASS\n")
    (blast_dir / "HA.blast.txt").write_text(
        "qHA\tACC_HA\tH5 influenza\t99.9\t1700\t1700\t1\t1700\t1\t1700\t0\t3000\n"
    )

    # NA: low confidence due to qcov < 90%
    (flags_dir / "NA.flag").write_text("PASS\n")
    (blast_dir / "NA.blast.txt").write_text(
        "qNA\tACC_NA\tN1 influenza\t99.9\t1200\t1350\t1\t1200\t1\t1200\t0\t2500\n"
    )

    # PB2: QC passed but no BLAST hit
    (flags_dir / "PB2.flag").write_text("PASS\n")
    (blast_dir / "PB2.blast.txt").write_text("")

    # Remaining segments: QC skipped
    for segment in ("PB1", "PA", "NP", "MP", "NS"):
        (flags_dir / f"{segment}.flag").write_text("FAIL\n")
        (blast_dir / f"{segment}.blast.txt").write_text("")

    output = tmp_path / "blast_top_hits.csv"
    fake = SimpleNamespace(
        wildcards=SimpleNamespace(sample="TEST_SAMPLE"),
        params=SimpleNamespace(min_identity=95.0, min_query_coverage=90.0),
        input=SimpleNamespace(
            blast_files=[str(blast_dir / f"{s}.blast.txt") for s in ("HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS")],
            flags=[str(flags_dir / f"{s}.flag") for s in ("HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS")],
        ),
        output=SimpleNamespace(csv=str(output)),
    )

    runpy.run_path(
        str(REPO_ROOT / "scripts/summarize_blast.py"),
        init_globals={"snakemake": fake},
    )

    with output.open(newline="", encoding="utf-8") as handle:
        rows = {row["segment"]: row for row in csv.DictReader(handle)}

    assert rows["HA"]["hit_status"] == "HIGH_CONFIDENCE"
    assert rows["NA"]["hit_status"] == "LOW_CONFIDENCE"
    assert rows["PB2"]["hit_status"] == "NO_HIT"
    assert rows["PB1"]["hit_status"] == "SKIPPED_QC"
    assert rows["HA"]["query_coverage"] == "100.000"
    assert rows["NA"]["query_coverage"] == "88.889"


def test_sample_summary_candidate_and_h5n1_review_flags(tmp_path):
    mod = load_module("sample_summary", "scripts/sample_summary.py")
    output = tmp_path / "summary.tsv"

    coverage_rows = [
        {
            "segment": "HA",
            "coverage_flag": "PASS",
            "contig": "A_HA_H5",
            "median_depth": "100",
            "candidate_count": "2",
        },
        {
            "segment": "NA",
            "coverage_flag": "FAIL",
            "contig": "A_NA_N1",
            "median_depth": "20",
            "candidate_count": "1",
        },
    ]

    fastplong = {
        "reads_before": 100,
        "bases_before": 1000,
        "q20_rate_before": 0.5,
        "q30_rate_before": 0.2,
        "reads_after": 90,
        "bases_after": 900,
        "q20_rate_after": 0.6,
        "q30_rate_after": 0.3,
        "reads_passed": 90,
        "reads_low_quality": 5,
        "reads_too_short": 5,
        "reads_with_adapters": 0,
    }

    mod.write_summary(
        output_path=output,
        sample="TEST_SAMPLE",
        fastplong=fastplong,
        coverage_rows=coverage_rows,
        blast_hits={"HA": "ACC_HA", "NA": "SKIPPED_QC"},
        h5n1_status="INDETERMINATE",
        genoflu_status="H5N1_INDETERMINATE",
        consensus_segments=1,
    )

    row = read_one_tsv(output)
    flags = set(row["review_flags"].split(";"))

    assert "h5n1_screen_indeterminate" in flags
    assert "coverage_failures" in flags
    assert "multiple_irma_candidates" in flags
    assert row["multiple_irma_candidate_segments"] == "HA"
    assert row["max_irma_candidate_count"] == "2"


def test_sample_summary_not_detected_is_not_itself_review_flag(tmp_path):
    mod = load_module("sample_summary_not_detected", "scripts/sample_summary.py")
    output = tmp_path / "summary.tsv"
    segments = ["HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS"]
    coverage_rows = [
        {
            "segment": segment,
            "coverage_flag": "PASS",
            "contig": "A_HA_H1" if segment == "HA" else "A_NA_N2" if segment == "NA" else segment,
            "median_depth": "100",
            "candidate_count": "1",
        }
        for segment in segments
    ]
    fastplong = {
        "reads_before": 100,
        "bases_before": 1000,
        "q20_rate_before": 0.5,
        "q30_rate_before": 0.2,
        "reads_after": 90,
        "bases_after": 900,
        "q20_rate_after": 0.6,
        "q30_rate_after": 0.3,
        "reads_passed": 90,
        "reads_low_quality": 5,
        "reads_too_short": 5,
        "reads_with_adapters": 0,
    }

    mod.write_summary(
        output_path=output,
        sample="TEST_SAMPLE",
        fastplong=fastplong,
        coverage_rows=coverage_rows,
        blast_hits={"HA": "ACC_HA", "NA": "ACC_NA"},
        h5n1_status="NOT_DETECTED",
        genoflu_status="H5N1_NOT_DETECTED",
        consensus_segments=8,
    )

    row = read_one_tsv(output)
    assert row["review_flags"] == "NONE"


def test_report_bundle_embeds_provenance(tmp_path):
    mod = load_module("build_report_bundle", "scripts/build_report_bundle.py")

    run_summary = tmp_path / "run_summary.html"
    run_summary.write_text("<html>run</html>", encoding="utf-8")
    sample_report = tmp_path / "S1.sample_summary.html"
    sample_report.write_text("<html>sample</html>", encoding="utf-8")
    provenance = tmp_path / "run_provenance.json"
    provenance.write_text(
        json.dumps({"workflow": {"name": "WINGS", "git_dirty": "false"}}),
        encoding="utf-8",
    )

    bundle = mod.build_bundle(run_summary, provenance, [sample_report])

    assert bundle["format"] == "WINGS_REPORT_BUNDLE"
    assert bundle["sample_count"] == 1
    assert bundle["provenance"]["filename"] == "run_provenance.json"
    assert bundle["provenance"]["data"]["workflow"]["name"] == "WINGS"


def test_sample_summary_preserves_metadata_and_protects_analytical_fields(tmp_path):
    mod = load_module("sample_summary_metadata", "scripts/sample_summary.py")
    output = tmp_path / "summary.tsv"

    segments = ["HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS"]
    coverage_rows = [
        {
            "segment": segment,
            "coverage_flag": "PASS",
            "contig": (
                "A_HA_H5" if segment == "HA"
                else "A_NA_N1" if segment == "NA"
                else segment
            ),
            "median_depth": "100",
            "candidate_count": "1",
        }
        for segment in segments
    ]

    fastplong = {
        "reads_before": 100,
        "bases_before": 1000,
        "q20_rate_before": 0.5,
        "q30_rate_before": 0.2,
        "reads_after": 90,
        "bases_after": 900,
        "q20_rate_after": 0.6,
        "q30_rate_after": 0.3,
        "reads_passed": 90,
        "reads_low_quality": 5,
        "reads_too_short": 5,
        "reads_with_adapters": 0,
    }

    metadata = {
        "sample_id": "TEST_SAMPLE",
        "host_common_name": "Mallard",
        "host_species": "Anas platyrhynchos",
        "collection_date": "2026-01-15",
        "state": "Arizona",
        "flyway": "Pacific",
        "custom_project_field": "sentinel-site-7",
        "review_flags": "metadata_must_not_override_analysis",
    }

    mod.write_summary(
        output_path=output,
        sample="TEST_SAMPLE",
        fastplong=fastplong,
        coverage_rows=coverage_rows,
        blast_hits={"HA": "ACC_HA", "NA": "ACC_NA"},
        h5n1_status="DETECTED",
        genoflu_status="COMPLETED",
        consensus_segments=8,
        metadata=metadata,
    )

    row = read_one_tsv(output)

    assert row["sample"] == "TEST_SAMPLE"
    assert row["sample_id"] == "TEST_SAMPLE"
    assert row["host_common_name"] == "Mallard"
    assert row["host_species"] == "Anas platyrhynchos"
    assert row["collection_date"] == "2026-01-15"
    assert row["state"] == "Arizona"
    assert row["flyway"] == "Pacific"
    assert row["custom_project_field"] == "sentinel-site-7"

    # Metadata cannot overwrite WINGS analytical fields.
    assert row["metadata_review_flags"] == "metadata_must_not_override_analysis"
    assert row["review_flags"] == "h5n1_screen_detected"
