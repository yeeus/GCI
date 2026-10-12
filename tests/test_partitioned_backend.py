"""Regression tests for opt-in read-name-partitioned BAM/PAF filtering.

Run with: python -m unittest discover -s tests -v
Requires pysam and the same dependencies as GCI.py.
"""
import gzip
import os
import shutil
import tempfile
import unittest

import numpy as np
import pysam

import GCI as gci


class PartitionedBackendTests(unittest.TestCase):
    @staticmethod
    def _make_reads(outdir):
        path = os.path.join(outdir, 'reads.bam')
        contigs = {'chrA': 800, 'chrB': 600}
        rows = []
        for i in range(80):
            chromosome = 'chrA' if i % 2 else 'chrB'
            start = (i * 17 + 11) % (contigs[chromosome] - 110)
            rows.append((f'read{i:03d}', chromosome, start, 100, 60 if i % 3 else 40))
        # Duplicate primary query names and reads crossing BAM split windows
        # exercise legacy ordering of per-window dictionary overwrites.
        rows.extend([
            ('dup', 'chrA', 50, 100, 60),
            ('dup', 'chrA', 200, 100, 40),
            ('chrB_only', 'chrB', 480, 100, 60),
        ])
        header = {
            'HD': {'VN': '1.6', 'SO': 'coordinate'},
            'SQ': [{'SN': chrom, 'LN': length} for chrom, length in contigs.items()],
        }
        with pysam.AlignmentFile(path, 'wb', header=header) as out:
            for name, chrom, start, size, mq in sorted(
                rows, key=lambda item: (list(contigs).index(item[1]), item[2])
            ):
                read = pysam.AlignedSegment()
                read.query_name = name
                read.query_sequence = 'A' * size
                read.flag = 0
                read.reference_id = list(contigs).index(chrom)
                read.reference_start = start
                read.mapping_quality = mq
                read.cigartuples = [(0, size)]
                read.set_tag('NM', 0)
                out.write(read)
        pysam.index(path)
        pafs = []
        for index in range(2):
            paf_path = os.path.join(outdir, f'reads{index}.paf')
            with open(paf_path, 'w') as out:
                for i in range(80):
                    if (i + index) % 3 == 0:
                        continue
                    chrom = 'chrA' if i % 2 else 'chrB'
                    start = (i * 17 + 11) % (contigs[chrom] - 110) + 5
                    out.write(
                        f'read{i:03d}\t100\t0\t100\t+\t{chrom}\t'
                        f'{contigs[chrom]}\t{start}\t{start + 100}\t'
                        f'100\t100\t{60 if i % 3 else 40}\n'
                    )
                out.write(
                    f'paf_only{index}\t100\t0\t100\t+\tchrA\t800\t'
                    f'100\t200\t100\t100\t60\n'
                )
            pafs.append(paf_path)
        return path, pafs

    def test_depths_match_legacy_with_zero_one_and_two_pafs(self):
        with tempfile.TemporaryDirectory() as outdir:
            bam, pafs = self._make_reads(outdir)
            for paf_count, threads, partitions in (
                (0, 1, 1), (1, 2, 5), (2, 4, 17),
            ):
                with self.subTest(pafs=paf_count, threads=threads, partitions=partitions):
                    legacy_prefix = f'old_{paf_count}'
                    new_prefix = f'new_{paf_count}'
                    legacy, legacy_targets = gci.filter(
                        paf_files=pafs[:paf_count], bam_files=[bam],
                        directory=outdir, prefix=legacy_prefix,
                        threads=threads, flank_len=15,
                    )
                    optimized, optimized_targets = gci.filter(
                        paf_files=pafs[:paf_count], bam_files=[bam],
                        directory=outdir, prefix=new_prefix,
                        threads=threads, flank_len=15,
                        filter_backend='partitioned',
                        scratch_dir=outdir,
                        partitions=partitions,
                    )
                    phase2c, phase2c_targets = gci.filter(
                        paf_files=pafs[:paf_count], bam_files=[bam],
                        directory=outdir, prefix=f'v2_{paf_count}',
                        threads=threads, flank_len=15,
                        filter_backend='partitioned-v2',
                        scratch_dir=outdir,
                        partitions=partitions,
                    )
                    self.assertEqual(legacy_targets, optimized_targets)
                    self.assertEqual(legacy_targets, phase2c_targets)
                    for chrom in legacy_targets:
                        np.testing.assert_array_equal(optimized[chrom], legacy[chrom])
                        np.testing.assert_array_equal(phase2c[chrom], legacy[chrom])
                    with gzip.open(os.path.join(outdir, legacy_prefix + '.depth.gz'), 'rb') as f:
                        a = f.read()
                    with gzip.open(os.path.join(outdir, new_prefix + '.depth.gz'), 'rb') as f:
                        b = f.read()
                    self.assertEqual(a, b)
                    with gzip.open(os.path.join(outdir, f'v2_{paf_count}.depth.gz'), 'rb') as f:
                        c = f.read()
                    self.assertEqual(a, c)
                    self.assertEqual(
                        gci.collapse_depth_range(legacy, -1, 0, 15),
                        gci.collapse_depth_range(optimized, -1, 0, 15),
                    )
                    self.assertEqual(
                        gci.collapse_depth_range(legacy, -1, 0, 15),
                        gci.collapse_depth_range(phase2c, -1, 0, 15),
                    )
            # Multiple BAM inputs use the same union/overlap semantics.
            second_bam = os.path.join(outdir, 'second.bam')
            shutil.copyfile(bam, second_bam)
            shutil.copyfile(bam + '.bai', second_bam + '.bai')
            legacy, _ = gci.filter(
                paf_files=pafs[:1], bam_files=[bam, second_bam],
                directory=outdir, prefix='legacy_two_bams',
                threads=2, flank_len=15,
            )
            optimized, _ = gci.filter(
                paf_files=pafs[:1], bam_files=[bam, second_bam],
                directory=outdir, prefix='partitioned_two_bams',
                threads=2, flank_len=15,
                filter_backend='partitioned',
                scratch_dir=outdir, partitions=7,
            )
            phase2c, _ = gci.filter(
                paf_files=pafs[:1], bam_files=[bam, second_bam],
                directory=outdir, prefix='phase2c_two_bams',
                threads=2, flank_len=15,
                filter_backend='partitioned-v2',
                scratch_dir=outdir, partitions=7,
            )
            for chrom in legacy:
                np.testing.assert_array_equal(legacy[chrom], optimized[chrom])
                np.testing.assert_array_equal(legacy[chrom], phase2c[chrom])
            self.assertFalse(any(name.startswith('gci-part-') for name in os.listdir(outdir)))


if __name__ == '__main__':
    unittest.main()
