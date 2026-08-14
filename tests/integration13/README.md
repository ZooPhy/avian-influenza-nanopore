# Phase 13 end-to-end release smoke test

Phase 13 is the highest-level compact WINGS integration test.

It starts from a raw Nanopore-style `*.fastq.gz` plus user metadata and ends at
the production `wings_report_bundle.wings` release artifact.

Production stages exercised directly:

1. metadata validation and per-sample extraction;
2. NanoPlot;
3. Porechop ABI;
4. fastplong;
5. seqtk read renaming;
6. IRMA normalization checkpoint;
7. coverage/QC;
8. Medaka model resolution;
9. real Medaka inference/consensus for the QC-passing HA segment;
10. BLAST and BLAST summary;
11. concatenated consensus;
12. H5N1 three-state screen;
13. GenoFLU gating;
14. production VADR rule/runtime boundary;
15. sample-summary TSV and HTML;
16. run-summary HTML;
17. run provenance;
18. final portable `.wings` report bundle.

Two boundaries are intentionally substituted:

- **IRMA assembly** is replaced by a deterministic IRMA-like project containing
  valid sorted/indexed BAMs and consensus FASTAs. The fixture is explicitly
  downstream of the real production `seqtk_rename` output.
- **VADR executable** is represented by a tiny local `v-annotate.pl` stub. The
  production VADR rule still executes its runtime selection and process handling.

GenoFLU itself is not stubbed. The fixture produces an `INDETERMINATE` H5N1
screen (HA passes, NA fails QC), so the production GenoFLU gate records
`H5N1_INDETERMINATE` without invoking GenoFLU.

The test also asserts that the production IRMA rule did not run by checking that
`irma/irma.log` was never created.

Run:

```bash
python -m pytest -q tests/integration13/test_end_to_end_release.py
```

Then:

```bash
python -m pytest -q
```
