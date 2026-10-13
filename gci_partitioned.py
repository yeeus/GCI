"""Disk-partitioned GCI BAM/PAF filtering backend.

Opt in via GCI.py --filter-backend partitioned. Retains legacy read-name
selection and coverage semantics without all read-level dictionaries in RAM.

Requires indexed coordinate-sorted BAMs, like the legacy backend.
Temporary files are cleaned on success and on Python exceptions.
"""

from collections import OrderedDict
from contextlib import nullcontext
from gci_depth_events import DepthEventAccumulator, INTERVAL_DTYPE
from multiprocessing import get_context
from pathlib import Path
import os
import tempfile
import time
import zlib

import numpy as np


class _BucketWriter:
    """Bound the number of simultaneously open temporary partition files."""

    def __init__(self, root, stem, buckets, max_open=256):
        self.root = Path(root)
        self.stem = stem
        self.buckets = buckets
        self.max_open = max_open
        self.handles = OrderedDict()

    def write(self, query, line):
        bucket = zlib.crc32(query.encode('utf-8')) % self.buckets
        if bucket in self.handles:
            fh = self.handles.pop(bucket)
        else:
            if len(self.handles) >= self.max_open:
                _, old = self.handles.popitem(last=False)
                old.close()
            fh = open(self.root / f'{self.stem}.{bucket:04d}.tsv', 'a')
        self.handles[bucket] = fh
        fh.write(line)

    def close(self):
        for fh in self.handles.values():
            fh.close()
        self.handles.clear()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


def _bam_task(task):
    """Read one legacy BAM window; write deduplicated mappings to a shard.

    Legacy overwrites multiple primary records of a query within a window,
    but ORs high-MAPQ membership. Maintain both.
    """
    import pysam

    task_id, target, left, right, bam_path, map_qual, clip_percent, iden_percent, mq_cutoff, root = task
    records = {}
    high_quality = set()
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for aln in bam.fetch(contig=target, start=left, stop=right, multiple_iterators=True):
            if not aln.is_mapped or aln.is_secondary or aln.is_supplementary or aln.mapping_quality < map_qual:
                continue
            cigar = aln.get_cigar_stats()[0]
            M, I, D, S, eq, X = (cigar[i] for i in (0, 1, 2, 4, 7, 8))
            mismatches = aln.get_tag('NM') - I - D
            if not (S / (M + eq + X + I + S) <= clip_percent):
                continue
            if not ((M + eq + X - mismatches) / (M + eq + X + I + D) >= iden_percent):
                continue
            query = aln.query_name
            records[query] = (aln.reference_name, aln.reference_start, aln.reference_end, aln.query_length)
            if aln.mapping_quality >= mq_cutoff:
                high_quality.add(query)

    output = Path(root) / f'bam-task-{task_id:06d}.tsv'
    with open(output, 'w') as out:
        for query, (chrom, start, end, qlen) in records.items():
            out.write(f'{query}\t{chrom}\t{start}\t{end}\t{qlen}\t{int(query in high_quality)}\n')
    return str(output), len(records)


def _merged_alignment_properties(alns, x, y):
    # Identical to legacy merge_alns_properties(), including overlaps.
    ordered = sorted([a[x:y + 1] for a in alns])
    lo, hi = ordered[0]
    intervals = []
    length = 0
    for left, right in ordered:
        if hi >= left:
            if hi < right:
                hi = right
        else:
            intervals.append((hi - lo, lo, hi))
            length += hi - lo
            lo, hi = left, right
    intervals.append((hi - lo, lo, hi))
    length += hi - lo
    intervals.sort(key=lambda item: item[0], reverse=True)
    return length, intervals[0][1], intervals[0][2]


def _paf_best_target(synteny):
    result = {}
    for query, by_target in synteny.items():
        by_score = {}
        for target, alns in by_target.items():
            aligned, _, _ = _merged_alignment_properties(alns, 1, 2)
            score = (sum(a[-1] for a in alns) / len(alns)) * (aligned / alns[0][0])
            _, start, end = _merged_alignment_properties(alns, 3, 4)
            by_score[target] = (score, start, end, alns[0][0])
        best_target = sorted(by_score, key=lambda t: (by_score[t][0], t), reverse=True)[0]
        score, start, end, qlen = by_score[best_target]
        result[query] = (best_target, start, end, qlen)
    return result


def _partition_worker(task):
    """Resolve one read-name bucket against all input mapping files."""
    bucket, root, n_paf, n_bam, targets, map_qual, mq_cutoff, identity_min, ovlp_min, binary_intervals, shard_plan = task
    root = Path(root)
    high_quality = set()
    synteny = {}
    paf_maps = [{} for _ in range(n_paf)]
    # Sharded inputs are read in original file / byte-range order.
    # The Phase 2B/C path continues using one bucket file with file IDs.
    if shard_plan is None:
        paf_inputs = [(None, root / f'paf.{bucket:04d}.tsv', True)]
    else:
        from gci_partition_fast import paf_bucket_fragments
        paf_inputs = [(file_id, path, False)
                      for file_id, path in paf_bucket_fragments(root, bucket, shard_plan[0])]
    current_file = 0
    for fixed_file_id, paf_path, has_file_id in paf_inputs:
        if not paf_path.exists():
            continue
        with open(paf_path) as src:
            for line in src:
                parts = line.rstrip('\n').split('\t')
                file_id = int(parts[0]) if has_file_id else fixed_file_id
                while current_file < file_id:
                    paf_maps[current_file] = _paf_best_target(synteny)
                    current_file += 1
                paf = parts[1:] if has_file_id else parts
                target = paf[5]
                if target not in targets:
                    continue
                qname = paf[0]
                qlen, qs, qe = int(paf[1]), int(paf[2]), int(paf[3])
                ts, te = int(paf[7]), int(paf[8])
                identity = int(paf[9]) / int(paf[10])
                mapq = int(paf[11])
                if mapq >= map_qual and identity >= identity_min:
                    synteny.setdefault(qname, {}).setdefault(target, []).append(
                        (qlen, qs, qe, ts, te, identity)
                    )
                    if mapq >= mq_cutoff:
                        high_quality.add(qname)
    # Also snapshots files without any records in this bucket.
    for file_id in range(current_file, n_paf):
        paf_maps[file_id] = _paf_best_target(synteny)

    del synteny
    bam_maps = []
    bam_fragments = [[] for _ in range(n_bam)]
    if shard_plan is None:
        for file_id in range(n_bam):
            bam_fragments[file_id].append(root / f'bam{file_id}.{bucket:04d}.tsv')
    else:
        from gci_partition_fast import bam_bucket_fragments
        for file_id, path in bam_bucket_fragments(root, bucket, shard_plan[1]):
            bam_fragments[file_id].append(path)
    for file_id in range(n_bam):
        reads = {}
        for bam_path in bam_fragments[file_id]:
            if not bam_path.exists():
                continue
            with open(bam_path) as src:
                for line in src:
                    qname, chrom, start, end, qlen, high = line.rstrip('\n').split('\t')
                    reads[qname] = (chrom, int(start), int(end), int(qlen))
                    if high == '1':
                        high_quality.add(qname)
        bam_maps.append(reads)

    files = paf_maps + bam_maps
    if len(files) > 1:
        common = set.intersection(*(set(f) for f in files))
        final = high_quality | common
        selected = {q: v for q, v in files[0].items() if q in final}
        for mapping in files[1:]:
            for q, segment in mapping.items():
                if q in selected:
                    previous = selected[q]
                    if segment[0] == previous[0]:
                        overlap = min(segment[2], previous[2]) - max(segment[1], previous[1])
                        if overlap / segment[-1] < ovlp_min:
                            del selected[q]
                        else:
                            selected[q] = (previous[0], max(segment[1], previous[1]), min(segment[2], previous[2]))
                    else:
                        del selected[q]
                elif q in high_quality:
                    selected[q] = (segment[0], segment[1], segment[2])
    else:
        selected = files[0]

    if binary_intervals:
        out_path = root / f'intervals.{bucket:04d}.bin'
        with open(out_path, 'wb') as out:
            chunk = []
            for segment in selected.values():
                chunk.append((targets[segment[0]], segment[1], segment[2]))
                if len(chunk) >= 65536:
                    np.asarray(chunk, dtype=INTERVAL_DTYPE).tofile(out)
                    chunk.clear()
            if chunk:
                np.asarray(chunk, dtype=INTERVAL_DTYPE).tofile(out)
    else:
        out_path = root / f'intervals.{bucket:04d}.tsv'
        with open(out_path, 'w') as out:
            for segment in selected.values():
                out.write(f'{segment[0]}\t{segment[1]}\t{segment[2]}\n')
    return str(out_path), len(selected)


def partitioned_filter(paf_files, bam_files, prefix, map_qual, mq_cutoff,
                       iden_percent, clip_percent, ovlp_percent, flank_len,
                       directory, force, log_reads_type, chrs_list, threads,
                       scratch_dir, partitions, write_depth, accumulation='direct',
                       input_mode='serial', input_shards=16):
    """Produce legacy-compatible (depths, target_lengths) with bounded partitions."""
    import pysam

    output_path = Path(directory) / f'{prefix}.depth.gz'
    if output_path.exists() and not force:
        raise SystemExit(f'ERROR!!! The file "{output_path}" exists. Please use -f')
    if threads < 1 or partitions < 1:
        raise ValueError('threads and partitions must be positive')
    if accumulation not in ('direct', 'events'):
        raise ValueError(f'Unknown accumulation mode: {accumulation}')
    if input_mode not in ('serial', 'parallel'):
        raise ValueError(f'Unknown input partition mode: {input_mode}')
    if input_shards < 1:
        raise ValueError('input_shards must be positive')
    if input_mode == 'parallel' and accumulation != 'events':
        raise ValueError('Parallel input partitioning requires event accumulation')
    binary_intervals = accumulation == 'events'
    print(f'Filtering {log_reads_type} alignment files ...', flush=True)

    with pysam.AlignmentFile(bam_files[0], 'rb') as bam:
        targets = dict(zip(bam.references, bam.lengths))
    if chrs_list:
        targets = {t: n for t, n in targets.items() if t in chrs_list}
    depths = {t: np.zeros(n, dtype=np.uint32) for t, n in targets.items()}
    t0 = time.perf_counter()
    profile = os.environ.get('GCI_PROFILE') == '1'

    with tempfile.TemporaryDirectory(prefix='gci-part-', dir=scratch_dir) as root:
        if input_mode == 'parallel':
            from gci_partition_fast import partition_paf_parallel
            paf_group_counts = partition_paf_parallel(
                paf_files, root, partitions, threads, input_shards
            )
        else:
            paf_group_counts = None
            # Existing Phase 2B/C PAF demultiplexing is unchanged.
            with _BucketWriter(root, 'paf', partitions) as writer:
                for i, path in enumerate(paf_files):
                    with open(path) as source:
                        for line in source:
                            query = line.split('\t', 1)[0]
                            writer.write(query, str(i) + '\t' + line)
        if profile:
            print(f'[GCI_PARTITION] type={log_reads_type} stage=paf_partition '
                  f'wall={time.perf_counter()-t0:.3f}s', flush=True)
        t_bam = time.perf_counter()

        if input_mode == 'parallel':
            from gci_partition_fast import partition_bam_parallel
            bam_group_counts = partition_bam_parallel(
                bam_files, targets, root, partitions, threads, input_shards,
                map_qual, clip_percent, iden_percent, mq_cutoff
            )
        else:
            bam_group_counts = None
            # Existing Phase 2B/C ordered BAM window redistribution is unchanged.
            for file_id, bam_path in enumerate(bam_files):
                tasks = []
                for target, length in targets.items():
                    step = 1 + (length - 1) // threads
                    for left in range(0, length, step):
                        tasks.append((len(tasks), target, left, min(left + step, length),
                                      bam_path, map_qual, clip_percent,
                                      iden_percent, mq_cutoff, root))
                with _BucketWriter(root, f'bam{file_id}', partitions) as writer:
                    if tasks:
                        with get_context('fork').Pool(processes=min(threads, len(tasks))) as pool:
                            for path, count in pool.imap(_bam_task, tasks, chunksize=1):
                                with open(path) as source:
                                    for line in source:
                                        writer.write(line.split('\t', 1)[0], line)
                                os.remove(path)
        if profile:
            print(f'[GCI_PARTITION] type={log_reads_type} stage=bam_partition '
                  f'wall={time.perf_counter()-t_bam:.3f}s', flush=True)
        t_merge = time.perf_counter()
        target_ids = {target: index for index, target in enumerate(targets)}
        worker_tasks = [(b, root, len(paf_files), len(bam_files),
                         target_ids if binary_intervals else set(targets),
                         map_qual, mq_cutoff, iden_percent, ovlp_percent,
                         binary_intervals,
                         (paf_group_counts, bam_group_counts)
                         if input_mode == 'parallel' else None)
                        for b in range(partitions)]
        accepted = 0
        accumulator_ctx = (
            DepthEventAccumulator(root, targets, flank_len=flank_len)
            if binary_intervals else nullcontext()
        )
        with accumulator_ctx as accumulator:
            with get_context('fork').Pool(processes=min(threads, partitions)) as pool:
                for path, count in pool.imap_unordered(_partition_worker, worker_tasks, chunksize=1):
                    accepted += count
                    if binary_intervals:
                        accumulator.add_file(path)
                    else:
                        with open(path) as source:
                            for line in source:
                                chrom, start, end = line.rstrip('\n').split('\t')
                                start = int(start) + flank_len
                                end = int(end) - flank_len
                                depths[chrom][start:end + 1] += 1
                    os.remove(path)
            if profile:
                print(f'[GCI_PARTITION] type={log_reads_type} '
                      f'stage={"filter_and_stage" if binary_intervals else "filter_and_accumulate"} '
                      f'wall={time.perf_counter()-t_merge:.3f}s accepted={accepted}', flush=True)
            if binary_intervals:
                t_prefix = time.perf_counter()
                accumulator.materialize(depths)
                if profile:
                    print(f'[GCI_PARTITION] type={log_reads_type} stage=depth_prefix_sum '
                          f'wall={time.perf_counter()-t_prefix:.3f}s '
                          f'intervals={accumulator.interval_count} '
                          f'events={accumulator.event_count}', flush=True)

    print(f'Filtering {log_reads_type} alignment files done!!!', flush=True)
    print(f'Writing depths into "{directory}/{prefix}.depth.gz" ...', flush=True)
    write_depth(directory, prefix, depths, threads)
    print('Writing depths done!!!\n\n', flush=True)
    return depths, targets
