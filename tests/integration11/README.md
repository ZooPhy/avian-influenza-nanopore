# Phase 11 production run-level summary integration

Phase 11 validates aggregation and reporting across multiple WINGS samples using
the real production `run_summary_html` rule and `scripts/run_summary.qmd`.

This phase deliberately starts at the run-reporting boundary. Earlier phases
already validate production normalization, coverage/QC, Medaka, BLAST, and
sample-summary generation. Here, deterministic per-sample outputs are supplied
for four samples representing every genome-completeness category:

| Sample | Passing segments | Expected genome status |
|---|---:|---|
| `complete` | 8 | Complete |
| `near_complete` | 7 | Near-complete |
| `partial` | 3 | Partial |
| `failed` | 0 | Failed |

The fixtures also span all three H5N1 screening states and include explicit
Medaka stage status files so the run-level Medaka aggregation is tested.

The production run-summary rule must generate:

- `run_summary/run_summary.html`
- `run_summary/run_summary.tsv`
- `run_summary/samples_requiring_review.tsv`

Assertions cover genome-status derivation, review selection, H/NA subtype
extraction from contig names, H5N1 reporting, Medaka status aggregation,
metadata preservation, and rendered HTML content.

Run:

```bash
python -m pytest -q tests/integration11/test_run_summary.py
```

Then:

```bash
python -m pytest -q
```
