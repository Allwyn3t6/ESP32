"""Offline Pan-Tompkins-style R-peak detector and scorer (ANSI/AAMI EC57-style 150 ms match window)."""
from __future__ import annotations

import numpy as np
from scipy.signal import butter, find_peaks, sosfiltfilt


def detect_r_peaks(x: np.ndarray, fs: float) -> np.ndarray:
    """Return R-peak times (s), sub-sample refined by parabolic interpolation."""
    sos = butter(2, [5, 15], btype="band", fs=fs, output="sos")
    bp = sosfiltfilt(sos, x)
    der = np.gradient(bp)
    w = int(0.15 * fs)
    mwi = np.convolve(der ** 2, np.ones(w) / w, mode="same")
    seg = int(10 * fs)
    thr = np.empty_like(mwi)
    for s in range(0, len(mwi), seg):
        thr[s:s + seg] = 0.3 * np.percentile(mwi[max(0, s - seg):s + 2 * seg], 99.5)
    pk, _ = find_peaks(mwi, distance=int(0.25 * fs))
    pk = pk[mwi[pk] > thr[pk]]
    hp = sosfiltfilt(butter(2, 0.5, btype="high", fs=fs, output="sos"), x)
    half = int(0.075 * fs)
    out = []
    for p in pk:
        a, b = max(1, p - half), min(len(x) - 1, p + half)
        i = a + int(np.argmax(hp[a:b]))
        if 1 <= i < len(x) - 1:
            y0, y1, y2 = hp[i - 1], hp[i], hp[i + 1]
            den = y0 - 2 * y1 + y2
            d = 0.5 * (y0 - y2) / den if den != 0 else 0.0
            out.append((i + np.clip(d, -0.5, 0.5)) / fs)
    return np.unique(np.round(np.array(out), 6))


def score(true_t: np.ndarray, det_t: np.ndarray, tol: float = 0.15):
    """Greedy nearest matching. Returns TP, FP, FN, timing errors (s) of TPs."""
    used = np.zeros(len(det_t), bool)
    errs = []
    if len(det_t) == 0:
        return 0, 0, len(true_t), np.array([])
    for t in true_t:
        i = int(np.searchsorted(det_t, t))
        cands = [j for j in (i - 1, i) if 0 <= j < len(det_t) and not used[j]]
        if not cands:
            continue
        j = min(cands, key=lambda j: abs(det_t[j] - t))
        if abs(det_t[j] - t) <= tol:
            used[j] = True
            errs.append(det_t[j] - t)
    tp = len(errs)
    return tp, int((~used).sum()), len(true_t) - tp, np.array(errs)
