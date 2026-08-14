# Phase 7 checkpoint-aware IRMA normalization integration

This phase exercises Snakemake's checkpoint reevaluation around the WINGS IRMA
normalization interface without running full IRMA.

The fixture creates an IRMA-like project containing:

- two complete HA candidates with different coverage medians;
- one NA FASTA-only candidate;
- one PB2 BAM-only candidate;
- five missing segments.

The test runs the production `scripts/normalize_irma_outputs.py` inside a real
Snakemake `checkpoint`, then uses `checkpoint.get()` in downstream input
functions to discover only normalized segments that are truly `READY`.

Expected normalization states:

| Segment | Expected |
|---|---|
| HA | READY, MULTIPLE_CANDIDATES; higher-depth candidate selected |
| NA | FASTA_ONLY |
| PB2 | BAM_ONLY |
| PB1/PA/NP/MP/NS | MISSING |

The downstream dynamically discovered file set must contain only HA.

This tests the structural behavior that the production WINGS Snakefile relies
on while avoiding the expensive IRMA assembly step.

Run:

```bash
python -m pytest -q tests/integration7/test_checkpoint_irma.py
```

Then:

```bash
python -m pytest -q
```
