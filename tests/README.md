# WINGS regression tests

This first regression layer tests stable workflow behavior without running IRMA,
Medaka, BLAST databases, containers, or Quarto. Tests use small synthetic inputs
and monkeypatch expensive external computation where needed.

Covered behaviors:

- segment QC: upper-length WARNING remains PASS when hard QC passes
- segment QC: low breadth fails
- multiple IRMA candidates remain a review condition rather than a hard QC failure
- BLAST summary states: HIGH_CONFIDENCE, LOW_CONFIDENCE, NO_HIT, SKIPPED_QC
- sample-summary review flags for H5N1 INDETERMINATE and candidate ambiguity
- H5N1 NOT_DETECTED is not itself a review flag
- portable `.wings` bundle embeds run-level provenance

Run from the repository root:

```bash
python -m pytest -q tests/test_regression.py
```

If pytest is not installed in the active Snakemake environment:

```bash
mamba install -n snakemake_env -c conda-forge pytest
```

These tests intentionally do not invoke the full Snakemake DAG. A later
integration layer should add a tiny prepared normalized-IRMA fixture and a CI
workflow.
