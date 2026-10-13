# GCI v1.1 — Release notes

**Date:** 2026-10-13  
**Previous release:** [v1.0](https://github.com/yeeus/GCI/releases/tag/v1.0)

GCI v1.1 promotes the validated Phase 2C filtering engine to the **only default implementation** and consolidates the maintenance updates made since v1.0. The GCI score interpretation and published output layout remain unchanged.

## Main changes

### 1. Faster, lower-memory filtering and depth processing

- **Read-name partitioning:** BAM and PAF records are partitioned by deterministic read-name hashing, allowing PAF best-target resolution and multi-mapper consistency filtering in bounded per-partition dictionaries rather than a monolithic whole-genome read dictionary.
- **Alignment filtering semantics preserved:** MAPQ and identity cutoffs, high-quality query retention, best-target scoring, PAF alignment block merging, cross-mapper target agreement and overlap thresholds retain their established behavior. Multiple PAF files retain the existing cumulative best-target convention, and BAM-window overwrite order is preserved.
- **Compact binary intervals:** Accepted alignment intervals use chromosome IDs and fixed-width start/end coordinates instead of large text intermediate files.
- **Difference events and prefix sums:** Coverage is calculated with signed interval-boundary events, then materialized in bounded 8-Mb windows into the original `uint32` depth arrays. Prefix-sum carry is preserved across window boundaries.
- **More efficient depth writing:** Fixed-size multiprocessing worker pools, faster repeated-depth text generation and ordered concatenation of gzip members remove much of the previous depth-writing bottleneck.
- **Lazy combined HiFi + ONT depth:** The combined track is generated as an on-demand maximum of HiFi and ONT depths rather than a fully materialized third genome-wide array. Each source track is gap-masked before combination.
- **Automatically managed temporary files:** Read partitions and depth events are written under a temporary `gci-part-*` subdirectory of the configured output directory, then removed when processing finishes or on Python exceptions.
- **Opt-in diagnostics:** Detailed `[GCI_PARTITION]` stage timers are now off by default; enable them with `GCI_PROFILE=1`.

### 2. Simplified CLI

- **Removed `--filter-backend`.** Phase 2C filtering is always active; the original filtering implementation and alternative development backends are not shipped in the v1.1 main source.
- **Removed `--scratch-dir`.** `-d PATH` specifies both the final output directory and parent directory for temporary partition files.
- **Retained `--partitions INT` with default 256.** Advanced users can tune the number of deterministic read-name partitions.
- **Removed `--input-shards`.** It belonged to the experimental Phase 2D backend.
- **`-t / --threads` remains available** for worker count and BAM-window partitioning; changes in `-t` can affect ambiguous repeated-query edge cases and should be validated when used for a new dataset.
- **Version string** updated to `GCI version 1.1`.

### 3. Main-repository improvements incorporated since v1.0

- Improved compatibility with BAM alignments using `--eqx` CIGAR conventions (mismatch/equality operations).
- Fixed an alignment-related `AttributeError` edge case and addressed depth-writer file-descriptor exhaustion.
- Improved plotting tools and plot appearance, including work on low-depth overlays, background handling, and efficient plotting of large datasets.
- Added and refined optional regional plotting/depth-mean functionality and related argument explanations.
- Updated workflow documentation for trio-based and Hi-C–phased assemblies, mapper inputs, MAPQ interpretation and stand-alone plotting.
- Added example MH63 results and supplementary plotting materials, with documentation and formatting cleanups.

These items reflect repository changes already present on the main branch after the v1.0 tag; the optional utility tools themselves are not rewritten by the Phase 2C filtering upgrade.

## Installation

For direct-script use, put **all three runtime Python files in the same directory**:

1. [GCI.py](https://github.com/yeeus/GCI/blob/main/GCI.py)
2. [gci_partitioned.py](https://github.com/yeeus/GCI/blob/main/gci_partitioned.py)
3. [gci_depth_events.py](https://github.com/yeeus/GCI/blob/main/gci_depth_events.py)

Dependencies include Python 3, `pysam`, `numpy`, `biopython`, and `matplotlib`. Indexed, coordinate-sorted BAM inputs are required. The existing Bioconda build may not contain the v1.1 code until repackaged.

```bash
python GCI.py -r assembly.fa \
  --hifi hifi.winnowmap.bam hifi.minimap2.paf \
  --nano ont.winnowmap.bam ont.minimap2.paf \
  -t 16 --partitions 256 -d output -o sample -f
```

No backend selection or scratch-directory flag is needed. Ensure `-d` has enough free space and adequate small-file I/O performance.

## MF2 full-genome validation

Complete MF2 maternal-haplotype assembly, mapped HiFi and ONT reads (each: winnowmap BAM + minimap2 PAF), with 16 workers:

| Metric | Earlier baseline | Phase 2C backend |
| --- | ---: | ---: |
| Elapsed wall time | 1,525.83 s (25m25.83s) | **470.08 s (7m50.08s)** |
| Maximum RSS | ~35.3 GiB | **~23.0 GiB** |
| Relative throughput | baseline | **~3.25× faster** |

On that benchmark, GCI reports, issue BED files and the decompressed HiFi, ONT and combined depth files matched the earlier implementation. This is validation of the **Phase 2C computational backend**; the v1.1 cleanup and moving scratch beneath `-d` also merit a release smoke test. Results are workload- and filesystem-dependent and are not universal speed claims.

## Compatibility notes and known limits

- The published output file names, GCI calculation and filtering thresholds are retained.
- At least one indexed BAM is required for each supplied sequencing technology.
- `-d` now holds automatically cleaned temporary files alongside final output files; use suitable high-capacity storage. Shared storage metadata operations can affect runtime.
- The default 256 partitions is configurable, not necessarily optimal for every sequencing depth or reference size.
- Read-name best-target/tie-breaking behavior and BAM-window read precedence are intentionally unchanged rather than silently “corrected”.
- Timing diagnostics are available with `GCI_PROFILE=1`; normal output is unchanged and does not show `[GCI_PARTITION]` timers.

## Experimental changes not included

Phase 2D parallel PAF/BAM input sharding and its Hybrid variant were evaluated during development but did not consistently outperform Phase 2C on full MF2 data. Neither these backends nor the development tests and benchmarking scripts are included in the minimal v1.1 main runtime.
