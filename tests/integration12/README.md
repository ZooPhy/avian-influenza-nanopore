# Phase 12 production bundle / release-artifact integration

Phase 12 validates the final portable `.wings` release artifact.

Earlier phases validate the analytical and reporting stages. This phase
pre-seeds deterministic rendered run/sample HTML reports and then executes the
real production:

1. `run_provenance`
2. `wings_report_bundle`

The production bundle builder stores the run-summary HTML, every sample-summary
HTML document, provenance JSON, bundle format/version metadata, generation
timestamp, and sample count inside a single JSON `.wings` file.

Assertions cover:

- bundle format `WINGS_REPORT_BUNDLE`;
- bundle version `1`;
- expected two-sample manifest;
- complete embedded run/sample HTML payloads;
- production provenance embedded exactly as written to disk;
- WINGS workflow identity and sample count;
- QC and BLAST parameters;
- BLAST database manifest metadata;
- hashes for all recorded Conda environment YAML files;
- byte-for-byte preservation after copying the `.wings` file into a standalone
  release directory;
- successful parsing with only the copied `.wings` artifact present;
- recovery ("rehydration") of all embedded HTML reports without referring back
  to the original WINGS results tree.

Run:

```bash
python -m pytest -q tests/integration12/test_report_bundle.py
```

Then:

```bash
python -m pytest -q
```
