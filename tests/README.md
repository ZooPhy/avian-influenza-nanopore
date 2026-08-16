# WINGS regression tests

The portable regression test layer checks stable workflow behavior without
running IRMA, Medaka, BLAST databases, containers, or Quarto. Tests use small
synthetic inputs and monkeypatch expensive external computation where needed.

Covered behaviors include:

- segment QC: upper-length WARNING remains PASS when hard QC passes
- segment QC: low breadth fails
- multiple IRMA candidates remain a review condition rather than a hard QC failure
- BLAST summary states: HIGH_CONFIDENCE, LOW_CONFIDENCE, NO_HIT, SKIPPED_QC
- sample-summary review flags for H5N1 INDETERMINATE and candidate ambiguity
- H5N1 NOT_DETECTED is not itself a review flag
- missing influenza segments remain distinct from failed segments in sample summaries
- portable `.wings` bundles embed run-level provenance

Run the portable regression tests from the repository root:

```bash
python -m pytest -q \
  tests/test_regression.py \
  tests/test_sample_summary_segments.py
```

If pytest is not installed in the active Snakemake environment:

```bash
conda install -n wings_snakemake_new -c conda-forge pytest
```

These portable tests intentionally do not execute the complete Snakemake DAG or
external bioinformatics tools. Full workflow behavior is validated separately
with integration runs on supported platforms. The portable regression suite is
also suitable for continuous integration because it does not require large
reference databases, container execution, or sequencing datasets.
