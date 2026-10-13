"""Phase 2D PAF parallel sharding: boundary and order regressions.

This test does not require pysam.
Run: python -m unittest discover -s tests -v
"""
import os
import tempfile
import unittest
import zlib

from gci_partition_fast import partition_paf_parallel, paf_bucket_fragments


class Phase2DInputTests(unittest.TestCase):
    def test_parallel_paf_byte_ranges_preserve_read_record_order(self):
        with tempfile.TemporaryDirectory() as root:
            files = []
            original = []
            for file_id in range(2):
                path = os.path.join(root, f'input{file_id}.paf')
                rows = []
                for i in range(35000):
                    query = f'read{(i * 37 + file_id) % 199}'
                    rows.append(
                        f'{query}\t80\t0\t80\t+\tchr1\t200000\t'
                        f'{i % 1000}\t{i % 1000 + 80}\t80\t80\t60\n'
                    )
                if file_id:
                    rows[-1] = rows[-1].rstrip('\n')
                with open(path, 'w') as out:
                    out.write(''.join(rows))
                files.append(path)
                original.append(rows)

            bucket_count = 13
            groups = partition_paf_parallel(
                files, root, bucket_count, threads=4, input_shards=4
            )
            self.assertTrue(all(n > 1 for n in groups))

            for bucket in range(bucket_count):
                actual = []
                for file_id, path in paf_bucket_fragments(root, bucket, groups):
                    with open(path) as handle:
                        actual.extend((file_id, line) for line in handle)
                expected = [
                    (file_id, line if line.endswith('\n') else line + '\n')
                    for file_id, rows in enumerate(original)
                    for line in rows
                    if zlib.crc32(line.split('\t', 1)[0].encode('utf-8')) % bucket_count == bucket
                ]
                self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
