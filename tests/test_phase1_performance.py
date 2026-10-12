"""Small regression tests for the Phase 1 performance changes.

Run from the repository root with: python -m unittest discover -s tests -v
"""
import gzip
import os
import tempfile
import unittest
from unittest import mock

import numpy as np
import pysam

import GCI as gci


class Phase1PerformanceTests(unittest.TestCase):
    def test_depth_writer_preserves_uncompressed_format(self):
        depths = {
            'chrA': (np.arange(16401, dtype=np.uint32) % 101),
            'chrB': np.array([0, 1, 5, 100], dtype=np.uint32),
        }
        with tempfile.TemporaryDirectory() as outdir:
            gci.write_depth(outdir, 'regression', depths, threads=2)
            with gzip.open(os.path.join(outdir, 'regression.depth.gz'), 'rt') as handle:
                actual = handle.read()
        expected = ''.join(
            f'>{target}\n' + ''.join(f'{int(depth)}\n' for depth in values)
            for target, values in depths.items()
        )
        self.assertEqual(actual, expected)

    def test_two_type_merge_preserves_maximum_depth(self):
        hifi = {'chrA': np.array([0, 10, 7, 5], dtype=np.uint32)}
        nano = {'chrA': np.array([8, 2, 7, 11], dtype=np.uint32)}
        with tempfile.TemporaryDirectory() as outdir:
            with mock.patch.object(gci, 'write_depth') as writer:
                merged = gci.merge_two_type_depth(hifi, nano, directory=outdir)
        np.testing.assert_array_equal(merged['chrA'], [8, 10, 7, 11])
        self.assertEqual(merged['chrA'].dtype, np.dtype('uint32'))
        writer.assert_called_once()

    def test_bam_filter_depth_and_dtype(self):
        with tempfile.TemporaryDirectory() as outdir:
            bam_path = os.path.join(outdir, 'reads.bam')
            header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                      'SQ': [{'SN': 'chrA', 'LN': 200}]}
            with pysam.AlignmentFile(bam_path, 'wb', header=header) as handle:
                read = pysam.AlignedSegment()
                read.query_name = 'test_read'
                read.query_sequence = 'A' * 80
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 40
                read.mapping_quality = 60
                read.cigartuples = [(0, 80)]
                read.set_tag('NM', 0)
                handle.write(read)
            pysam.index(bam_path)

            depths, targets = gci.filter(
                bam_files=[bam_path], directory=outdir,
                prefix='bam_regression', threads=1, flank_len=15
            )
            self.assertEqual(targets, {'chrA': 200})
            self.assertEqual(depths['chrA'].dtype, np.dtype('uint32'))
            expected = np.zeros(200, dtype=np.uint32)
            expected[55:106] = 1  # Preserve the existing inclusive end convention.
            np.testing.assert_array_equal(depths['chrA'], expected)


if __name__ == '__main__':
    unittest.main()
