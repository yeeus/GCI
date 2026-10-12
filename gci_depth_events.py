"""Bounded-memory difference-event accumulation for the GCI partitioned backend.

Each accepted interval has a compact, fixed-width binary representation:
chromosome ID (signed int32), start and end (signed int64). Coordinates
represent the legacy untrimmed interval and are normalized with Python's
slice semantics before generating coverage events.

Event records are uint32: low bit 0 for +1 or 1 for -1, remaining bits store
the offset within a depth window. Window size must be <= 2**31 - 1.
"""

from collections import OrderedDict
from pathlib import Path
import numpy as np

INTERVAL_DTYPE = np.dtype([('chrom', '<i4'), ('start', '<i8'), ('end', '<i8')])


class DepthEventAccumulator:
    """Stage bounded binary events and materialize uint32 depths per window."""

    def __init__(self, directory, chromosomes, flank_len=15,
                 window_size=8_000_000, max_open=64):
        if not 1 <= window_size < (1 << 31):
            raise ValueError('window_size must lie in [1, 2**31 - 1]')
        self.directory = Path(directory)
        self.chromosomes = list(chromosomes.items())
        self.flank_len = flank_len
        self.window_size = window_size
        self.max_open = max_open
        self.handles = OrderedDict()
        self.interval_count = 0
        self.event_count = 0

    def _path(self, chrom_id, window):
        return self.directory / f'events.{chrom_id:04d}.{window:05d}.u32'

    def _handle(self, chrom_id, window):
        key = (chrom_id, window)
        if key in self.handles:
            return self.handles.pop(key)
        if len(self.handles) >= self.max_open:
            _, old = self.handles.popitem(last=False)
            old.close()
        return open(self._path(*key), 'ab')

    def _emit(self, chrom_id, positions, negative):
        if len(positions) == 0:
            return
        windows = positions // self.window_size
        offsets = (positions % self.window_size).astype(np.uint32)
        codes = (offsets << 1) | np.uint32(1 if negative else 0)
        for win in np.unique(windows):
            mask = windows == win
            stream = self._handle(chrom_id, int(win))
            self.handles[(chrom_id, int(win))] = stream
            codes[mask].tofile(stream)
            self.event_count += int(np.count_nonzero(mask))

    def add(self, records):
        """Add a batch of raw accepted intervals without storing it long-term."""
        if records.size == 0:
            return
        if records.dtype != INTERVAL_DTYPE:
            raise TypeError('Binary intervals must use INTERVAL_DTYPE')
        if np.any(records['chrom'] < 0) or np.any(records['chrom'] >= len(self.chromosomes)):
            raise ValueError('Interval refers to unknown chromosome ID')
        self.interval_count += len(records)
        for chrom_id in np.unique(records['chrom']):
            chrom_id = int(chrom_id)
            length = self.chromosomes[chrom_id][1]
            by_chr = records[records['chrom'] == chrom_id]
            # Reproduce depths[chrom][start:end + 1] += 1, including negative
            # indices and clipping at array boundaries.
            starts = by_chr['start'] + self.flank_len
            stops = by_chr['end'] - self.flank_len + 1
            starts = np.clip(np.where(starts < 0, starts + length, starts), 0, length)
            stops = np.clip(np.where(stops < 0, stops + length, stops), 0, length)
            valid = stops > starts
            if not np.any(valid):
                continue
            starts, stops = starts[valid], stops[valid]
            self._emit(chrom_id, starts, negative=False)
            # An end event exactly at the chromosome end would never be read.
            self._emit(chrom_id, stops[stops < length], negative=True)

    def add_file(self, path, batch_size=65536):
        path = Path(path)
        if path.stat().st_size % INTERVAL_DTYPE.itemsize:
            raise ValueError(f'Truncated binary interval file: {path}')
        with open(path, 'rb') as inp:
            while True:
                batch = np.fromfile(inp, dtype=INTERVAL_DTYPE, count=batch_size)
                if batch.size == 0:
                    break
                self.add(batch)

    def close(self):
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()

    def materialize(self, depths, event_batch_size=1_000_000):
        """Apply exact integer prefix sums, preserving cross-window coverage."""
        self.close()
        for chrom_id, (chrom, length) in enumerate(self.chromosomes):
            depth = depths[chrom]
            if len(depth) != length or depth.dtype != np.dtype('uint32'):
                raise ValueError('Depth track does not match its chromosome')
            carry = 0
            for base in range(0, length, self.window_size):
                width = min(self.window_size, length - base)
                path = self._path(chrom_id, base // self.window_size)
                if not path.exists():
                    if carry:
                        depth[base:base + width] = carry
                    continue
                diff = np.zeros(width, dtype=np.int64)
                with open(path, 'rb') as inp:
                    while True:
                        codes = np.fromfile(inp, dtype='<u4', count=event_batch_size)
                        if not codes.size:
                            break
                        offsets = (codes >> 1).astype(np.intp)
                        if np.any(offsets >= width):
                            raise ValueError(f'Out-of-window event in {path}')
                        deltas = 1 - (codes & 1).astype(np.int64) * 2
                        np.add.at(diff, offsets, deltas)
                balance = int(np.sum(diff, dtype=np.int64))
                np.cumsum(diff, dtype=np.int64, out=diff)
                if carry:
                    diff += carry
                depth[base:base + width] = diff
                carry += balance
                path.unlink()

    def __enter__(self):
        return self

    def __exit__(self, typ, val, tb):
        self.close()
