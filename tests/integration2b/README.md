# WINGS Phase 2B integration tests

Phase 2B exercises WINGS segment QC and BLAST summarization through Snakemake
without running IRMA, Medaka, or a real BLAST database.

The test harness uses the production scripts:

- `scripts/check_coverage.py`
- `scripts/coverage_table.py`
- `scripts/summarize_blast.py`

It creates tiny synthetic HA alignments and checks four important paths:

- `qc_pass`: complete 60x HA coverage -> QC `PASS` -> BLAST `HIGH_CONFIDENCE`
- `length_warning`: 1909-nt HA with complete 60x coverage and two IRMA candidates -> QC `PASS` with `WARNING` -> BLAST `LOW_CONFIDENCE` at 89.261% query coverage
- `low_breadth`: only 1500/1700 HA positions covered at 60x -> QC `FAIL` -> BLAST `SKIPPED_QC`
- `no_hit`: QC `PASS`, empty BLAST evidence -> `NO_HIT`

The BLAST step itself is mocked deliberately. The purpose here is to test the
Snakemake handoff from segment QC into the production BLAST-summary script, not
to rebuild or query the large influenza reference database.

Run from the repository root:

```bash
python -m pytest -q tests/integration2b/test_qc_blast_integration.py
```

Generated files are placed under `tests/integration2b/work/` and are ignored by Git.
