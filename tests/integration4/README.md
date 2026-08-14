# Phase 4 real BLAST smoke test

This test adds the first integration layer that executes real NCBI BLAST+.

It:

1. creates a deterministic synthetic 1700-nt HA reference and query;
2. builds a tiny nucleotide database with `makeblastdb`;
3. runs real `blastn` using the same outfmt and ranking logic as WINGS;
4. feeds the raw BLAST output into the production `scripts/summarize_blast.py`;
5. verifies an exact HA match becomes `HIGH_CONFIDENCE`;
6. verifies a QC-failing NA segment remains `SKIPPED_QC`.

The BLAST executable comes from the production `envs/blast.yaml` Conda environment.

Run:

```bash
python -m pytest -q tests/integration4/test_real_blast.py
```

Then:

```bash
python -m pytest -q
```
