"""Independent Phase 2C tests; no pysam needed."""
import os
import tempfile
import unittest

import numpy as np

from gci_depth_events import INTERVAL_DTYPE, DepthEventAccumulator


class DepthEventTests(unittest.TestCase):
    def compare(self, targets, rows, flank, window, batch_size=3):
        expected = {chrom: np.zeros(n, dtype=np.uint32) for chrom, n in targets.items()}
        for chrom_id, start, end in rows:
            chrom = list(targets)[chrom_id]
            expected[chrom][start + flank:end - flank + 1] += 1
        actual = {chrom: np.zeros(n, dtype=np.uint32) for chrom, n in targets.items()}
        with tempfile.TemporaryDirectory() as root:
            input_file = os.path.join(root, 'intervals.bin')
            np.array(rows, dtype=INTERVAL_DTYPE).tofile(input_file)
            with DepthEventAccumulator(root, targets, flank, window, max_open=2) as acc:
                acc.add_file(input_file, batch_size=batch_size)
                self.assertEqual(acc.interval_count, len(rows))
                acc.materialize(actual, event_batch_size=3)
            self.assertFalse(any(x.startswith('events.') for x in os.listdir(root)))
        for chrom in targets:
            np.testing.assert_array_equal(actual[chrom], expected[chrom])

    def test_cross_window_and_exact_slice_boundaries(self):
        rows = [
            (0, 0, 40), (0, 2, 2), (0, 11, 11),
            (0, -30, 15), (0, 10, 90), (0, 18, 25),
            (0, 19, 20), (0, 19, 21), (0, 0, 0),
            (1, 0, 40), (1, 5, 6), (1, -2, 18),
        ] + [(0, 0, 19)] * 100
        for flank in (0, 1, 4, 15):
            for window in (1, 2, 3, 5, 8, 19, 20, 21):
                with self.subTest(flank=flank, window=window):
                    self.compare({'chrA': 20, 'chrB': 31}, rows, flank, window)

    def test_randomized_against_legacy_uint32_slice(self):
        rng = np.random.default_rng(41)
        for _ in range(100):
            targets = {'chrA': int(rng.integers(1, 150)),
                       'chrB': int(rng.integers(1, 240))}
            rows = [(int(rng.integers(0, 2)), int(rng.integers(-500, 500)),
                     int(rng.integers(-500, 500)))
                    for _ in range(int(rng.integers(0, 240)))]
            self.compare(targets, rows, flank=int(rng.integers(0, 30)),
                         window=int(rng.integers(1, 90)))

    def test_invalid_binary_interval_rejected(self):
        with tempfile.TemporaryDirectory() as root:
            file = os.path.join(root, 'bad.bin')
            with open(file, 'wb') as out:
                out.write(b'bad')
            with DepthEventAccumulator(root, {'chrA': 10}) as accumulator:
                with self.assertRaisesRegex(ValueError, 'Truncated'):
                    accumulator.add_file(file)


if __name__ == '__main__':
    unittest.main()
