# Phase 9 production checkpoint -> Medaka -> BLAST integration

Phase 9 connects the major downstream WINGS stages in one production-style
chain using the real top-level `Snakefile`.

The fixture contains one deterministic 1700-nt HA segment at 60x nominal
coverage and a tiny BLAST database containing the exact synthetic reference.

The helper workflow only prepares:

- a valid IRMA-like project with sorted/indexed HA BAM and consensus FASTA;
- a real NCBI nucleotide BLAST database.

The production WINGS workflow then executes:

1. `normalize_irma_outputs` checkpoint;
2. `check_coverage`;
3. `resolve_medaka_model`;
4. real `medaka_inference`;
5. real `medaka_consensus`;
6. real `blastn`.

This intentionally targets `HA.blast.txt` rather than a whole-sample summary so
that only one segment enters neural-network inference.

Expected assertions include:

- HA normalization = READY / UNIQUE;
- QC = PASS with 60x median depth and breadth 1.000;
- FASTQ-derived Medaka selector is resolved;
- real `features.hdf` is non-empty;
- Medaka consensus status = SUCCESS / MEDAKA;
- BLAST output is non-empty;
- the top hit is the synthetic H5 reference at >=99% identity and >=99% query
  coverage.

Run:

```bash
python -m pytest -q tests/integration9/test_checkpoint_medaka_blast.py
```

Then run the complete suite:

```bash
python -m pytest -q
```
