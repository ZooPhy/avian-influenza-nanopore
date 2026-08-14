# WINGS Phase 2 Snakemake integration tests

This layer exercises **actual rules from the production `Snakefile`** without running IRMA, Medaka, BLAST, VADR, or Quarto.

The fixture generator creates three tiny synthetic cases:

- `it_detected`: QC-passing H5 + N1 -> `DETECTED`
- `it_not_detected`: informative non-H5 HA + N1 -> `NOT_DETECTED`
- `it_indeterminate`: failed HA QC -> `INDETERMINATE`, with a multiple-IRMA-candidate review condition

The integration test asks the production Snakefile to build each `*.sample_summary.tsv`. All other required inputs are pre-created. The H5N1 flag is intentionally absent, so Snakemake must wire:

```text
coverage_flags + coverage_stats
        |
        v
detect_h5n1
        |
        v
sample_summary
```

Only the production rules `detect_h5n1` and `sample_summary` are allowed in this test. This prevents accidental upstream IRMA/Medaka/BLAST execution.

Run from the repository root:

```bash
python -m pytest -q tests/integration/test_snakemake_integration.py
```

The generated workspace lives under `tests/integration/work/` and should not be committed.
