# Phase 5 Medaka smoke test

This phase intentionally does **not** run neural-network Medaka inference.

It validates the architecture-sensitive pieces around Medaka:

- the production `envs/medaka.yaml` can create an environment and exposes
  Medaka 2.2.2;
- the production `scripts/resolve_medaka_model.py` resolves the Oxford Nanopore
  `basecall_model_version_id` from FASTQ metadata;
- the resulting consensus and variant selectors are correct;
- a QC-failing segment produces the expected `SKIPPED_QC` state;
- a failed inference with fail-soft behavior preserves the QC-passing IRMA
  consensus as `IRMA_FALLBACK`.

This keeps Phase 5 deterministic and fast while testing the state machine that
protects the production workflow from Medaka failures.

Run:

```bash
python -m pytest -q tests/integration5/test_medaka_smoke.py
```

Then run the complete suite:

```bash
python -m pytest -q
```
