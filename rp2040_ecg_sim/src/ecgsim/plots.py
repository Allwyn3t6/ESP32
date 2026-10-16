"""Matplotlib figures (PNG, 200 dpi)."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COL = {"A": "#2a78d6", "B": "#eb6834", "C": "#1baf7a"}
LBL = {"A": "A: CPU polling + USB CDC", "B": "B: ADC DMA + USB CDC", "C": "C: ADC DMA + PIO/DMA link"}
plt.rcParams.update({"figure.dpi": 100, "savefig.dpi": 200, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": "#e4e3df",
                     "grid.linewidth": 0.6, "font.size": 9, "axes.titlesize": 10, "legend.frameon": False})


def _save(fig, path):
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def grouped_bars(agg, conds, metric, ylabel, title, path, log=False, scale=1.0):
    fig, ax = plt.subplots(figsize=(10, 4))
    w = 0.26
    x = np.arange(len(conds))
    for i, a in enumerate("ABC"):
        sub = agg[agg.arch == a].set_index("condition").reindex(conds)
        m = sub[f"{metric}_mean"].values * scale
        lo = sub[f"{metric}_ci_lo"].values * scale
        hi = sub[f"{metric}_ci_hi"].values * scale
        err = np.vstack([np.nan_to_num(m - lo), np.nan_to_num(hi - m)])
        mm = np.where(log, np.maximum(m, 1e-12), m) if log else m
        ax.bar(x + (i - 1) * w, mm, w - 0.02, color=COL[a], label=LBL[a], yerr=err, capsize=2,
               error_kw=dict(lw=0.8, ecolor="#52514e"))
    ax.set_xticks(x, [c.replace("_", "\n") for c in conds])
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if log:
        ax.set_yscale("log")
    ax.legend(ncol=3, loc="upper left", bbox_to_anchor=(0, -0.18))
    _save(fig, path)


def waveform_excerpt(raw, fs, path):
    fig, axs = plt.subplots(3, 1, figsize=(10, 6.5), sharex=True)
    t = np.arange(len(raw["A"]["ref"])) / fs
    sel = (t >= 2) & (t < 6)
    for ax, a in zip(axs, "ABC"):
        ax.plot(t[sel], raw[a]["ref"][sel], color="#8a8985", lw=1.2, label="noise-free reference")
        ax.plot(t[sel], raw[a]["v"][sel], color=COL[a], lw=0.9, label="received (reconstructed)")
        ax.set_ylabel("ADC input (V)")
        ax.set_title(LBL[a], loc="left")
        ax.legend(loc="upper right", ncol=2)
    axs[-1].set_xlabel("time (s)")
    _save(fig, path)


def interval_hist(raw, path):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
    e = raw["A"]["interval_err_us"]
    axs[0].hist(e, bins=200, color=COL["A"])
    axs[0].set_yscale("log")
    axs[0].set_xlabel("sampling-interval error (µs)")
    axs[0].set_ylabel("count")
    axs[0].set_title("A: CPU polling (trial 0, normal)")
    for a in "BC":
        axs[1].hist(raw[a]["interval_err_us"] * 1e3, bins=100, color=COL[a], alpha=0.6, label=LBL[a])
    axs[1].set_xlabel("sampling-interval error (ns)")
    axs[1].set_title("B and C: ADC CLKDIV pacing (identical)")
    axs[1].legend()
    _save(fig, path)


def latency_cdf(raws, path, key="lat_ms", title="End-to-end transport latency"):
    fig, axs = plt.subplots(1, len(raws), figsize=(4.2 * len(raws), 3.6), sharey=True)
    for ax, (cond, raw) in zip(np.atleast_1d(axs), raws.items()):
        for a in "ABC":
            x = np.sort(raw[a][key])
            ax.plot(x, 1 - np.arange(x.size) / x.size, color=COL[a], lw=2, label=LBL[a])
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("latency (ms)")
        ax.set_title(cond.replace("_", " "))
    np.atleast_1d(axs)[0].set_ylabel("P(latency > x)")
    np.atleast_1d(axs)[0].legend(loc="lower left")
    fig.suptitle(title)
    _save(fig, path)


def sensitivity(sdf, sweep, metrics, path, xlabel, logx=False):
    sub = sdf[sdf.sweep == sweep]
    fig, axs = plt.subplots(1, len(metrics), figsize=(4.0 * len(metrics), 3.4))
    for ax, (m, lab) in zip(np.atleast_1d(axs), metrics):
        for a in "ABC":
            g = sub[sub.arch == a].groupby("value")[m]
            mu, sd = g.mean(), g.std()
            n = g.count()
            h = 2.262 * sd / np.sqrt(n)   # t(0.975, 9)
            ax.errorbar(mu.index, mu.values, yerr=h.values, color=COL[a], marker="o", ms=4, lw=2,
                        capsize=2, label=LBL[a])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(lab)
        if logx:
            ax.set_xscale("symlog", linthresh=1e-9) if sub.value.min() == 0 else ax.set_xscale("log")
    np.atleast_1d(axs)[0].legend(loc="best")
    _save(fig, path)
