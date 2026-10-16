"""Host receiver: CRC check, sequence-number based reconstruction, gap/dup/reorder detection."""
from __future__ import annotations

import numpy as np

from . import packets as P


def reconstruct(items, n_total):
    """items: iterable of (deliver_time, bytes, ...) in delivery order.

    Returns dict with codes (float, NaN = missing), receive time per sample, counters."""
    codes = np.full(n_total, np.nan)
    rx_time = np.full(n_total, np.nan)
    got = np.zeros(n_total, bool)          # slot seen (either data or 'acquisition-missed' marker)
    acq_marker = np.zeros(n_total, bool)
    dup = ooo = crc_fail = 0
    max_seq = -1
    n_ok_packets = 0
    for it in items:
        t, b = it[0], it[1]
        parsed = P.parse_packet(b)
        if parsed is None:
            crc_fail += 1
            continue
        seq0, c = parsed
        n_ok_packets += 1
        idx = seq0 + np.arange(len(c))
        valid = idx < n_total
        idx, c = idx[valid], c[valid]
        seen = got[idx]
        dup += int(seen.sum())
        new = ~seen
        ooo += int((idx[new] < max_seq).sum())
        i2, c2 = idx[new], c[new]
        got[i2] = True
        miss = c2 == P.MISSING
        acq_marker[i2[miss]] = True
        codes[i2[~miss]] = c2[~miss]
        rx_time[i2] = t
        if idx.size:
            max_seq = max(max_seq, int(idx.max()))
    return dict(codes=codes, rx_time=rx_time, got=got, acq_marker=acq_marker,
                duplicates=dup, out_of_order=ooo, crc_fail=crc_fail, ok_packets=n_ok_packets)


def missing_runs(got: np.ndarray, N: int):
    """Number of lost packets inferred purely from sequence gaps (slot runs never seen)."""
    miss = ~got
    if not miss.any():
        return 0
    return int(np.ceil(miss.sum() / N))
