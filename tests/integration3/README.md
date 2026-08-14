# Phase 3 downstream smoke test

This integration test exercises the production `Snakefile` rules:

- `detect_h5n1`
- `sample_summary`

The fixture pre-seeds upstream results so the test does **not** run IRMA,
Medaka inference, BLAST, VADR, or Quarto.

It validates one coherent mixed-outcome sample in which HA passes as H5,
NA fails segment QC, the H5N1 screen becomes `INDETERMINATE`, review flags
are generated, and validated metadata is propagated into the final
`sample_summary.tsv`.

Run only Phase 3:

```bash
python -m pytest -q tests/integration3/test_downstream_smoke.py
```

Run the full suite:

```bash
python -m pytest -q
```
