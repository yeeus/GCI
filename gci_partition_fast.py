"""Parallel direct-to-bucket input partitioning for GCI Phase 2D.

Each PAF shard is a contiguous byte range of an input file. Each BAM shard
is a contiguous group of the original BAM scan windows. The consumer reads
the shards in their original order, preserving repeated read-name overwrite
semantics and cumulative multi-PAF best-target selection.
"""

from multiprocessing import get_context
from pathlib import Path
import os


def _partition_paf_chunk(task):
    """Split a byte interval by read name without splitting PAF lines."""
    file_id, path, group_id, start, stop, directory, partitions = task
    from gci_partitioned import _BucketWriter

    count = 0
    stem = f'paf{file_id}.group{group_id:03d}'
    with open(path, 'rb') as source, _BucketWriter(directory, stem, partitions) as writer:
        if start:
            source.seek(start - 1)
            if source.read(1) != b'\n':
                source.readline()
        else:
            source.seek(0)
        while source.tell() < stop:
            line = source.readline()
            if not line:
                break
            decoded = line.decode('utf-8')
            if not decoded.endswith('\n'):
                decoded += '\n'
            query = decoded.split('\t', 1)[0]
            writer.write(query, decoded)
            count += 1
    return count


def partition_paf_parallel(paths, directory, partitions, threads, input_shards):
    """Return PAF group counts, preserving input-file and byte-range order."""
    counts = []
    tasks = []
    for i, path in enumerate(paths):
        length = os.path.getsize(path)
        groups = min(input_shards, max(1, (length + 1048575) // 1048576))
        counts.append(groups)
        for group in range(groups):
            tasks.append((i, path, group, length * group // groups,
                          length * (group + 1) // groups,
                          directory, partitions))
    if tasks:
        with get_context('fork').Pool(processes=min(threads, len(tasks))) as pool:
            list(pool.imap(_partition_paf_chunk, tasks, chunksize=1))
    return counts


def _partition_bam_group(task):
    """Filter BAM windows and emit directly to read-name buckets."""
    bam_id, group_id, tasks, directory, partitions = task
    from gci_partitioned import _BucketWriter
    import pysam

    stem = f'bam{bam_id}.group{group_id:03d}'
    accepted = 0
    with _BucketWriter(directory, stem, partitions) as writer:
        for (target, left, right, path, map_qual, clip_percent,
             iden_percent, mq_cutoff) in tasks:
            records = {}
            high_quality = set()
            with pysam.AlignmentFile(path, 'rb') as bam:
                for aln in bam.fetch(contig=target, start=left, stop=right,
                                     multiple_iterators=True):
                    if (not aln.is_mapped or aln.is_secondary or
                            aln.is_supplementary or aln.mapping_quality < map_qual):
                        continue
                    cigar = aln.get_cigar_stats()[0]
                    M, I, D, S, eq, X = (cigar[i] for i in (0, 1, 2, 4, 7, 8))
                    mismatches = aln.get_tag('NM') - I - D
                    if not S / (M + eq + X + I + S) <= clip_percent:
                        continue
                    if not (M + eq + X - mismatches) / (M + eq + X + I + D) >= iden_percent:
                        continue
                    query = aln.query_name
                    records[query] = (aln.reference_name, aln.reference_start,
                                      aln.reference_end, aln.query_length)
                    if aln.mapping_quality >= mq_cutoff:
                        high_quality.add(query)
            for query, (chrom, start, end, qlen) in records.items():
                writer.write(query,
                             f'{query}\t{chrom}\t{start}\t{end}\t{qlen}\t'
                             f'{int(query in high_quality)}\n')
            accepted += len(records)
    return accepted


def partition_bam_parallel(bam_files, targets, directory, partitions,
                           threads, input_shards, map_qual, clip_percent,
                           iden_percent, mq_cutoff):
    """Return BAM group counts, avoiding the main-process redistribution."""
    counts = []
    for bam_id, bam_path in enumerate(bam_files):
        tasks = []
        for target, length in targets.items():
            step = 1 + (length - 1) // threads
            for left in range(0, length, step):
                tasks.append((target, left, min(left + step, length), bam_path,
                              map_qual, clip_percent, iden_percent, mq_cutoff))
        groups = min(input_shards, len(tasks)) if tasks else 0
        counts.append(groups)
        if not groups:
            continue
        chunks = []
        for group_id in range(groups):
            left = len(tasks) * group_id // groups
            right = len(tasks) * (group_id + 1) // groups
            chunks.append((bam_id, group_id, tasks[left:right], directory, partitions))
        with get_context('fork').Pool(processes=min(threads, groups)) as pool:
            list(pool.imap(_partition_bam_group, chunks, chunksize=1))
    return counts


def paf_bucket_fragments(root, bucket, counts):
    root = Path(root)
    for file_id, groups in enumerate(counts):
        for group in range(groups):
            path = root / f'paf{file_id}.group{group:03d}.{bucket:04d}.tsv'
            if path.exists():
                yield file_id, path


def bam_bucket_fragments(root, bucket, counts):
    root = Path(root)
    for bam_id, groups in enumerate(counts):
        for group in range(groups):
            path = root / f'bam{bam_id}.group{group:03d}.{bucket:04d}.tsv'
            if path.exists():
                yield bam_id, path
