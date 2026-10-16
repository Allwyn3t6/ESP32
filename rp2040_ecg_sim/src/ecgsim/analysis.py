"""Aggregation, confidence intervals and paired statistical tests."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

METRICS = [
    "interval_err_rms_us", "interval_err_std_us", "interval_err_maxabs_us", "abs_time_err_rms_us",
    "acq_sample_loss_rate", "rx_missing_rate", "rx_transport_missing_rate", "duplicates", "out_of_order",
    "adc_fifo_overflow", "device_buffer_overflow_pkts", "receiver_buffer_overflow_pkts",
    "packet_loss_rate", "corrupted_packets", "crc_detected", "corruption_detection_rate",
    "lat_mean_ms", "lat_median_ms", "lat_p95_ms", "lat_p99_ms", "lat_max_ms", "lat_jitter_std_ms",
    "link_lat_mean_ms", "link_lat_p99_ms", "link_lat_jitter_std_ms",
    "rmse_mv", "rmse_vs_afe_mv", "rmse_filled_mv", "snr_db",
    "r_sensitivity", "r_ppv", "r_timing_err_rms_ms", "r_timing_err_maxabs_ms", "cpu_util_est",
]

# metrics used for paired hypothesis tests (lower is better unless listed in HIGHER_BETTER)
TEST_METRICS = ["interval_err_rms_us", "acq_sample_loss_rate", "rx_missing_rate", "packet_loss_rate",
                "device_buffer_overflow_pkts", "receiver_buffer_overflow_pkts",
                "lat_mean_ms", "lat_p99_ms", "lat_jitter_std_ms", "link_lat_jitter_std_ms",
                "rmse_vs_afe_mv", "rmse_filled_mv", "r_sensitivity", "r_timing_err_rms_ms"]
HIGHER_BETTER = {"r_sensitivity", "r_ppv", "snr_db"}


def ci95(x):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    if x.size < 2:
        return np.nan, np.nan
    h = stats.t.ppf(0.975, x.size - 1) * x.std(ddof=1) / np.sqrt(x.size)
    return x.mean() - h, x.mean() + h


def aggregate(df: pd.DataFrame, by=("condition", "arch")) -> pd.DataFrame:
    out = []
    for key, g in df.groupby(list(by), sort=False):
        row = dict(zip(by, key if isinstance(key, tuple) else (key,)))
        row["n_trials"] = len(g)
        for m in METRICS:
            if m not in g:
                continue
            x = g[m].astype(float)
            lo, hi = ci95(x)
            row[f"{m}_mean"] = x.mean()
            row[f"{m}_std"] = x.std(ddof=1)
            row[f"{m}_ci_lo"], row[f"{m}_ci_hi"] = lo, hi
        out.append(row)
    return pd.DataFrame(out)


def holm(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    order = np.argsort(pv)
    m = pv.size
    adj = np.empty(m)
    run = 0.0
    for r, i in enumerate(order):
        run = max(run, (m - r) * pv[i])
        adj[i] = min(1.0, run)
    out[ok] = adj
    return out


def paired_tests(df: pd.DataFrame, pairs=(("A", "B"), ("B", "C"), ("A", "C"))) -> pd.DataFrame:
    """Paired (same trial seed => same ECG and acquisition draws) comparisons.
    Wilcoxon signed-rank test (two-sided, zero differences dropped) + paired-t CI of the mean difference.
    Holm-Bonferroni correction over all tests in the table."""
    rows = []
    for cond, g in df.groupby("condition", sort=False):
        piv = {a: ga.set_index("trial") for a, ga in g.groupby("arch")}
        for x, y in pairs:
            for m in TEST_METRICS:
                a, b = piv[x][m].astype(float), piv[y][m].astype(float)
                d = (a - b).dropna()
                r = dict(condition=cond, comparison=f"{x} vs {y}", metric=m, n=len(d),
                         mean_x=a.mean(), mean_y=b.mean(), mean_diff=d.mean())
                lo, hi = ci95(d) if len(d) > 1 else (np.nan, np.nan)
                r["diff_ci_lo"], r["diff_ci_hi"] = lo, hi
                if len(d) == 0 or np.allclose(d, 0, atol=1e-12, rtol=0):
                    r["p_wilcoxon"] = 1.0
                    r["effect_rank_biserial"] = 0.0
                else:
                    res = stats.wilcoxon(d, zero_method="wilcox", alternative="two-sided")
                    r["p_wilcoxon"] = float(res.pvalue)
                    nz = d[d != 0]
                    rk = stats.rankdata(np.abs(nz))
                    r["effect_rank_biserial"] = float((rk[nz > 0].sum() - rk[nz < 0].sum()) / rk.sum())
                lower_better = m not in HIGHER_BETTER
                r["better"] = ("tie" if r["p_wilcoxon"] >= 1.0 else
                               (x if (r["mean_diff"] < 0) == lower_better else y))
                rows.append(r)
    out = pd.DataFrame(rows)
    out["p_holm"] = holm(out["p_wilcoxon"].values)
    out["significant_0.05"] = out["p_holm"] < 0.05
    out.loc[~out["significant_0.05"], "better"] = "no significant difference"
    return out
