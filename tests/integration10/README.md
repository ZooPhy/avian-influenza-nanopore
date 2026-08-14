# Phase 10 production sample-summary integration

Phase 10 validates the reporting boundary using the real production WINGS
`Snakefile`, `scripts/sample_summary.py`, and `scripts/sample_summary.qmd`.

The fixture creates eight normalized-ready IRMA-like segments. Only HA passes
the production 50x coverage gate; the other seven segments are held at 20x.
That keeps real Medaka neural-network inference limited to HA while still
allowing the production coverage table to represent all eight segments.

The production chain exercised is:

1. `normalize_irma_outputs`
2. `check_coverage`
3. `coverage_table`
4. `resolve_medaka_model`
5. `medaka_inference`
6. `medaka_consensus`
7. `blastn`
8. `summarize_blast`
9. `concat_consensus`
10. `detect_h5n1`
11. `sample_summary`
12. `sample_summary_html`

The test deliberately pre-seeds deterministic inputs that are not the reporting
focus: fastplong JSON, per-sample metadata, GenoFLU status, VADR log, and Medaka
VCF/status files.

Key reporting assertions include:

- arbitrary metadata is propagated;
- a metadata field named `segments_pass` becomes `metadata_segments_pass`;
- analytical `segments_pass` remains authoritative;
- HA PASS and seven failed segments are summarized correctly;
- H5N1 is `INDETERMINATE` because NA fails segment QC;
- review flags include incomplete segment QC and coverage failures;
- the real HA BLAST hit reaches the summary;
- the Quarto HTML report renders and contains metadata, H5N1 status, and BLAST
  content.

Run:

```bash
python -m pytest -q tests/integration10/test_sample_summary.py
```

Then:

```bash
python -m pytest -q
```
