"""Small regression tests for the Phase 1 performance changes.

Run from the repository root with: python -m unittest discover -s tests -v
"""
import gzip
import io
import os
from contextlib import redirect_stdout
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

    def test_depth_writer_handles_long_runs_and_frequent_changes(self):
        depths = {
            'chrA': np.array(
                [0] * 17000 + [5] * 17000 + [4294967295] * 17,
                dtype=np.uint32
            ),
            'chrB': np.tile(np.array([0, 1], dtype=np.uint32), 17000),
        }
        with tempfile.TemporaryDirectory() as outdir:
            gci.write_depth(outdir, 'runs', depths, threads=2)
            with gzip.open(os.path.join(outdir, 'runs.depth.gz'), 'rt') as handle:
                actual = handle.read()
        expected = ''.join(
            f'>{target}\n' + ''.join(f'{int(value)}\n' for value in values)
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

    def test_lazy_combined_depth_matches_eager_and_gap_masks(self):
        hifi = {
            'chrA': np.array([9, 2, 0, 7, 0, 0, 4, 6], dtype=np.uint32),
            'chrB': np.arange(33001, dtype=np.uint32) % 7,
        }
        nano = {
            'chrA': np.array([1, 8, 5, 3, 0, 1, 2, 7], dtype=np.uint32),
            'chrB': np.arange(33001, dtype=np.uint32) % 11,
        }
        gaps = {'chrA': [(1, 3)], 'chrB': [(16382, 16390)]}
        gci.merge_gaps_depths(hifi, gaps)
        gci.merge_gaps_depths(nano, gaps)
        with tempfile.TemporaryDirectory() as outdir:
            lazy = gci.merge_two_type_depth(
                hifi, nano, prefix='lazy', directory=outdir,
                threads=4, lazy=True
            )
            self.assertNotIsInstance(lazy['chrA'], np.ndarray)
            for target in hifi:
                expected = np.maximum(hifi[target], nano[target])
                np.testing.assert_array_equal(lazy[target][:], expected)
                np.testing.assert_array_equal(lazy[target][2:20], expected[2:20])
                self.assertEqual(len(lazy[target]), len(expected))
            with gzip.open(os.path.join(outdir, 'lazy.depth.gz'), 'rt') as f:
                observed = f.read()
        expected_text = ''.join(
            f'>{target}\n' +
            ''.join(f'{int(value)}\n' for value in np.maximum(hifi[target], nano[target]))
            for target in hifi
        )
        self.assertEqual(observed, expected_text)
        expected_intervals = gci.collapse_depth_range(
            {target: np.maximum(hifi[target], nano[target]) for target in hifi},
            -1, 0, 0,
        )
        self.assertEqual(
            gci.collapse_depth_range(lazy, -1, 0, 0),
            expected_intervals,
        )

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

            # A matching PAF from another mapper should restrict the
            # accepted interval without changing the filtering rule.
            paf_path = os.path.join(outdir, 'reads.paf')
            with open(paf_path, 'w') as handle:
                handle.write('test_read\t80\t0\t80\t+\tchrA\t200\t45\t125\t80\t80\t60\n')
            combined, _ = gci.filter(
                paf_files=[paf_path], bam_files=[bam_path],
                directory=outdir, prefix='bam_paf_regression',
                threads=1, flank_len=15
            )
            expected_combined = np.zeros(200, dtype=np.uint32)
            expected_combined[60:106] = 1
            np.testing.assert_array_equal(combined['chrA'], expected_combined)

            with mock.patch.dict(os.environ, {'GCI_PROFILE': '1'}):
                capture = io.StringIO()
                with redirect_stdout(capture):
                    profiled, _ = gci.filter(
                        paf_files=[paf_path], bam_files=[bam_path],
                        directory=outdir, prefix='profile_regression',
                        threads=1, flank_len=15,
                    )
            np.testing.assert_array_equal(profiled['chrA'], combined['chrA'])
            profile_log = capture.getvalue()
            for stage in (
                'depth_init', 'paf_parse', 'paf_best_target',
                'bam_pool_map', 'bam_results_merge',
                'cross_mapper_filter', 'depth_accumulation', 'depth_write',
            ):
                self.assertIn(f'stage={stage}', profile_log)


    def test_issue_intervals_preserve_original_boundaries(self):
        depth = np.full(100, 3, dtype=np.uint32)
        depth[15:25] = 0
        depth[30:40] = 0
        depth[75:85] = 0
        # The original implementation suppresses runs ending within
        # the first flank_len positions of the scanned interval.
        observed = gci.collapse_depth_range({'chrA': depth}, -1, 0, 15)
        self.assertEqual(observed['chrA'], [(30, 40), (75, 85)])

        regional = np.array([0, 0, 4, 0, 0, 5, 0], dtype=np.uint32)
        observed = gci.collapse_depth_range(
            {'chrA': regional}, -1, 0, 0, 10
        )
        self.assertEqual(observed['chrA'], [(10, 12), (13, 15), (16, 17)])

    def test_issue_interval_spanning_chunk_boundary(self):
        depth = np.full(8_000_100, 2, dtype=np.uint32)
        depth[7_999_990:8_000_040] = 0
        observed = gci.collapse_depth_range({'chrA': depth}, -1, 0, 15)
        self.assertEqual(observed['chrA'], [(7_999_990, 8_000_040)])


if __name__ == '__main__':
    unittest.main()
