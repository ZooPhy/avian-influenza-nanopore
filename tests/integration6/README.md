# Phase 6 real Medaka inference smoke test

This is the first WINGS integration test that executes real neural-network
Medaka inference.

The fixture is intentionally tiny:

- one deterministic synthetic 1400-nt HA-like contig;
- 40 synthetic Oxford Nanopore-style reads derived from that contig;
- approximately 1.2% substitution noise;
- a FASTQ header carrying
  `basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.0.0`.

The test then:

1. runs the production `scripts/resolve_medaka_model.py`;
2. aligns the reads with `minimap2` and creates a sorted/indexed BAM with
   `samtools`;
3. runs real `medaka inference`;
4. requires a non-empty real `features.hdf`;
5. runs real `medaka sequence`;
6. requires a non-empty polished consensus and `SUCCESS` status records.

The commands intentionally mirror the production WINGS Medaka commands, while
avoiding the IRMA checkpoint dependency that would otherwise make this a full
pipeline run.

The production Medaka environment pins `medaka==2.2.2`.

Run:

```bash
python -m pytest -q tests/integration6/test_real_medaka.py
```

The first run can take longer if a Medaka model must be resolved or downloaded.

Then run the complete suite:

```bash
python -m pytest -q
```
