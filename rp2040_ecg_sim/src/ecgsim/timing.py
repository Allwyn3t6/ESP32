"""Busy-interval / stall-window primitives used by the CPU and link models."""
from __future__ import annotations

import numpy as np


def lognormal_mean(rng, mean, sigma, size):
    """Lognormal samples parameterised by their arithmetic mean."""
    mu = np.log(mean) - sigma ** 2 / 2
    return rng.lognormal(mu, sigma, size)


def lognormal_median(rng, median, sigma, size):
    return rng.lognormal(np.log(median), sigma, size)


class Intervals:
    """A set of disjoint, sorted [start, end) intervals with O(log n) queries."""

    def __init__(self, start=None, end=None):
        self.start = np.asarray([] if start is None else start, dtype=float)
        self.end = np.asarray([] if end is None else end, dtype=float)

    @classmethod
    def fifo_queue(cls, arrivals, durations):
        """Jobs served FIFO by one server (Lindley recursion, vectorised):
        end_i = S_i + max_{j<=i}(a_j - S_{j-1}). Returns non-overlapping busy intervals."""
        a = np.asarray(arrivals, float)
        d = np.asarray(durations, float)
        if a.size == 0:
            return cls()
        S = np.cumsum(d)
        end = S + np.maximum.accumulate(a - (S - d))
        return cls(end - d, end)

    @classmethod
    def poisson(cls, rng, rate, T, dur_fn, queue=True):
        """Poisson-arriving windows of random duration over [0, T]."""
        if rate <= 0 or T <= 0:
            return cls()
        n = rng.poisson(rate * T)
        arr = np.sort(rng.uniform(0, T, n))
        dur = dur_fn(n)
        if queue:
            return cls.fifo_queue(arr, dur).merged()
        return cls(arr, arr + dur).merged()

    def merged(self):
        """Merge touching/overlapping intervals so a single push suffices."""
        if self.start.size == 0:
            return self
        o = np.argsort(self.start)
        s, e = self.start[o], self.end[o]
        e_run = np.maximum.accumulate(e)
        new = np.ones(s.size, bool)
        new[1:] = s[1:] > e_run[:-1]
        gid = np.cumsum(new) - 1
        ns = s[new]
        ne = np.zeros(ns.size)
        np.maximum.at(ne, gid, e_run)
        return Intervals(ns, ne)

    def end_if_inside(self, t):
        """If t lies inside an interval return that interval's end, else t."""
        t = np.asarray(t, float)
        if self.start.size == 0:
            return t.copy() if t.ndim else float(t)
        i = np.searchsorted(self.start, t, side="right") - 1
        ic = np.clip(i, 0, None)
        inside = (i >= 0) & (t < self.end[ic])
        return np.where(inside, self.end[ic], t)

    def push(self, t: float) -> float:
        """Scalar version of end_if_inside (merged intervals -> one lookup)."""
        if self.start.size == 0:
            return t
        i = int(np.searchsorted(self.start, t, side="right")) - 1
        if i >= 0 and t < self.end[i]:
            return float(self.end[i])
        return t

    @property
    def total(self):
        return float(np.sum(self.end - self.start))
