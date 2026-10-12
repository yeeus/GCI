# GCI disk-partitioned filtering (experimental)

This backend is **opt-in** and does not replace the verified legacy backend.

## Run

From the directory containing GCI.py:

```bash
python GCI.py -r /path/assembly.fa \
  --hifi /path/hifi.bam /path/hifi.paf \
  --nano /path/ont.bam /path/ont.paf \
  --filter-backend partitioned \
  --scratch-dir /path/to/node-local/scratch \
  --partitions 256 -t 16 -d output -o sample -f
```

- `--filter-backend legacy` is the default. The new backend must be requested explicitly.
- `--scratch-dir` must exist and have sufficient free disk space. Prefer node-local SSD/NVMe scratch, not a shared filesystem. When omitted, Python uses the system temporary directory (`TMPDIR` or `/tmp`).
- `--partitions` controls hash buckets for read-name grouping; 256 is the default. More buckets reduce peak per-worker read dictionaries but increase partition file count and filesystem overhead.
- `-t` controls the BAM processing and read-partition worker counts. Start with 16 on real data. Try 32 or 64 only after output equivalence is verified.
- Set `GCI_PROFILE=1` to print `[GCI_PARTITION]` timing for PAF partitioning, BAM partitioning, and partition filtering/depth accumulation.

## Algorithm

1. Hash each PAF record by read name into temporary bucket files.
2. Filter BAM with the existing chromosome-window scan, but serialize each window's query dictionary into a temporary file rather than returning large dictionaries over multiprocessing pipes. Redistribute those records to the same read-name bucket, preserving the legacy window overwrite order.
3. Process independent buckets in parallel. The worker reproduces the legacy PAF scoring, cumulative handling of multiple PAF files, high-MAPQ union, cross-mapper target/overlap checks, and final read selection.
4. Write the resulting intervals into temporary files; in the parent, reuse the exact legacy `depths[target][start:end+1] += 1` accumulation with flank trimming, followed by the verified depth writer.

The temporary working directory is removed automatically on normal exit or Python exception.

## Validation

Before running on MF2, run all tests:

```bash
python -m unittest discover -s tests -v
```

Compare the new and previous outputs, including decompressed depth data:

```bash
diff old/sample.gci new/sample.gci
for type in hifi nano two_type; do
  diff "old/sample_\${type}.0.depth.bed" "new/sample_\${type}.0.depth.bed"
  cmp <(gzip -cd "old/sample_\${type}.depth.gz") <(gzip -cd "new/sample_\${type}.depth.gz")
done
```

## Limits and expected behavior

- First release optimizes **read-level memory and parallel cross-mapper processing**, not chromosome-wide coverage storage.
- All input reads are partitioned on disk; this incurs substantial extra I/O and temporary disk usage. End-to-end speedup is not guaranteed, especially on networked scratch.
- Initial PAF demultiplexing and BAM shard redistribution are still sequential; their costs must be measured before more optimization.
- Currently Linux only for multiprocessing (`fork`).
- Coordinate-sorted, indexed BAM files are required, as in legacy.
- The current implementation may still use substantial RAM during `uint32` genome-wide depth accumulation and when processing unusually skewed read-name buckets.
- The legacy code has several edge-case semantics (PAF cumulative across multiple files, BAM window conflict precedence, inclusive depth slice end) that are deliberately preserved here. These are **not** silently corrected.
- Results have been checked in local synthetic tests and a new regression test has been added; **full-data MF2 runtime/memory and output equality remain to be verified on your cluster**.

Do not merge into `main` until the regressions and a complete MF2 comparison have passed.
