# Phase 8 production checkpoint -> coverage/QC integration

This phase runs the **production WINGS Snakefile** from the normalized IRMA
checkpoint through segment QC and the final coverage table.

Full IRMA assembly is intentionally not rerun. A helper Snakemake file first
builds a deterministic IRMA-like project containing valid sorted/indexed BAMs
and consensus FASTAs for seven recovered segments. **MP is deliberately absent**
so the production workflow must represent a genuinely unrecovered segment as
`MISSING` rather than failing during DAG construction. The production workflow
then executes only:

- `normalize_irma_outputs` (real checkpoint)
- `check_coverage` (real production rule/script)
- `coverage_table` (real production rule/script)

The fixture exercises seven recovered segments plus one genuinely absent segment:

| Segment | Fixture | Expected |
|---|---|---|
| HA | 1700 nt, 60x | PASS |
| NA | 1400 nt, 20x | coverage FAIL |
| PB2 | 2450 nt, 60x | length WARNING, overall PASS |
| PB1 | 2100 nt, 60x | length FAIL |
| PA | 2200 nt, 50 Ns, 60x | N-content FAIL |
| NP | 1500 nt, 60x | PASS |
| MP | absent from IRMA-like project | MISSING |
| NS | 850 nt, 60x | PASS |

Coverage is computed by the production code with
`samtools depth -aa -q 0 -Q 0`, using the WINGS thresholds of 50x median
depth and 0.95 breadth at that depth.

The test additionally verifies that the normalized MP FASTA/BAM do not exist,
that the manifest records MP as `MISSING`, that QC emits an MP `MISSING` flag,
and that the final coverage table still contains all eight influenza segments.

Run:

```bash
python -m pytest -q tests/integration8/test_checkpoint_coverage.py
```

Then:

```bash
python -m pytest -q
```
