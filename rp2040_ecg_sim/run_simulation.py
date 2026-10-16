#!/usr/bin/env python3
"""Main entry point: runs the full Monte Carlo study, sensitivity sweeps, statistics and plots.

Usage:
    python run_simulation.py                 # full study (config.json)
    python run_simulation.py --quick         # 3 trials x 60 s smoke run (results_quick/)
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import platform
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
from ecgsim import analysis, plots  # noqa: E402
from ecgsim.experiment import apply_overrides, run_trial  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def _job(args):
    cfg, label, extra, trial, keep = args
    rows, raw = run_trial(cfg, trial, keep_raw=keep)
    for r in rows:
        r.update(label)
        r.update(extra)
    return rows, (label, raw) if keep else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(HERE, "config.json"))
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--skip-sensitivity", action="store_true")
    args = ap.parse_args()
    base = json.load(open(args.config))
    out = os.path.join(HERE, "results_quick" if args.quick else "results")
    os.makedirs(os.path.join(out, "csv"), exist_ok=True)
    os.makedirs(os.path.join(out, "plots"), exist_ok=True)
    n_trials = 3 if args.quick else base["n_trials"]
    if args.quick:
        base["duration_s"] = 60
    t0 = time.time()

    # ---------------- main experiment
    jobs = []
    for cond, ov in base["conditions"].items():
        cfg = apply_overrides(base, ov)
        for tr in range(n_trials):
            jobs.append((cfg, {"condition": cond}, {}, tr, tr == 0))
    with Pool(base.get("n_workers", os.cpu_count())) as pool:
        res = pool.map(_job, jobs, chunksize=1)
    rows = [r for rr, _ in res for r in rr]
    raws = {lab["condition"]: raw for _, x in res if x for lab, raw in [x]}
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, "csv", "trial_metrics.csv"), index=False)
    agg = analysis.aggregate(df)
    agg.to_csv(os.path.join(out, "csv", "aggregate_metrics.csv"), index=False)
    tests = analysis.paired_tests(df)
    tests.to_csv(os.path.join(out, "csv", "paired_tests.csv"), index=False)
    print(f"main experiment: {len(df)} rows in {time.time() - t0:.0f} s")

    # ---------------- sensitivity sweeps
    sdf = None
    if not args.skip_sensitivity:
        sc = base["sensitivity"]
        sjobs = []
        for name, sw in sc.items():
            if not isinstance(sw, dict):
                continue
            cfg0 = apply_overrides(base, base["conditions"][sw["base"]])
            cfg0 = apply_overrides(cfg0, sw["fixed"])
            cfg0["duration_s"] = 30 if args.quick else sc["duration_s"]
            for v in sw["values"]:
                ov = {sw["param"]: v}
                if sw["fixed"].get("_tie_host_stall_to_param"):
                    ov["host.stall_median_s"] = v
                cfg = apply_overrides(cfg0, ov)
                for tr in range(2 if args.quick else sc["n_trials"]):
                    sjobs.append((cfg, {"sweep": name}, {"value": v, "param": sw["param"]}, tr, False))
        with Pool(base.get("n_workers", os.cpu_count())) as pool:
            sres = pool.map(_job, sjobs, chunksize=1)
        sdf = pd.DataFrame([r for rr, _ in sres for r in rr])
        sdf.to_csv(os.path.join(out, "csv", "sensitivity_trials.csv"), index=False)
        analysis.aggregate(sdf, by=("sweep", "value", "arch")).to_csv(
            os.path.join(out, "csv", "sensitivity_aggregate.csv"), index=False)
        print(f"sensitivity: {len(sdf)} rows, total {time.time() - t0:.0f} s")

    # ---------------- plots
    P = os.path.join(out, "plots")
    conds = list(base["conditions"])
    gb = plots.grouped_bars
    gb(agg, conds, "interval_err_rms_us", "RMS interval error (µs, log)", "Acquisition sampling-interval error (RMS)",
       f"{P}/01_interval_error_rms.png", log=True)
    gb(agg, conds, "interval_err_maxabs_us", "max |interval error| (µs, log)", "Acquisition sampling-interval error (max)",
       f"{P}/02_interval_error_max.png", log=True)
    gb(agg, conds, "acq_sample_loss_rate", "samples lost at acquisition (%)", "Acquisition sample loss",
       f"{P}/03_acq_sample_loss.png", scale=100)
    gb(agg, conds, "rx_missing_rate", "missing at receiver (%)", "Receiver missing-sample rate (all causes)",
       f"{P}/04_rx_missing.png", scale=100)
    gb(agg, conds, "device_buffer_overflow_pkts", "packets per 10-min record", "Device TX/DMA buffer overflows",
       f"{P}/05_device_buffer_overflow.png")
    gb(agg, conds, "receiver_buffer_overflow_pkts", "packets per 10-min record", "Receiver buffer overflows (no flow control)",
       f"{P}/06_receiver_buffer_overflow.png")
    gb(agg, conds, "lat_p99_ms", "p99 latency (ms, log)", "End-to-end transport latency, 99th percentile",
       f"{P}/07_latency_p99.png", log=True)
    gb(agg, conds, "link_lat_jitter_std_ms", "std of link latency (ms, log)", "Device-to-receiver link latency jitter",
       f"{P}/08_link_latency_jitter.png", log=True)
    gb(agg, conds, "rmse_vs_afe_mv", "RMSE vs AFE output (mV)", "Waveform error vs ideally-sampled AFE output (received samples)",
       f"{P}/09_rmse_vs_afe.png")
    gb(agg, conds, "rmse_filled_mv", "RMSE incl. gap interpolation (mV)", "Waveform reconstruction error vs noise-free reference",
       f"{P}/10_rmse_filled.png")
    gb(agg, conds, "r_timing_err_rms_ms", "RMS R-peak timing error (ms)", "R-peak timing error",
       f"{P}/11_rpeak_timing.png")
    gb(agg, conds, "r_sensitivity", "sensitivity", "R-peak detection sensitivity", f"{P}/12_rpeak_sensitivity.png")
    if "normal" in raws:
        plots.waveform_excerpt(raws["normal"], base["fs"], f"{P}/00_waveform_excerpt.png")
        plots.interval_hist(raws["normal"], f"{P}/13_interval_error_hist.png")
        sel = {c: raws[c] for c in ("normal", "bursty_delays", "receiver_delay") if c in raws}
        plots.latency_cdf(sel, f"{P}/14_latency_ccdf.png")
        plots.latency_cdf(sel, f"{P}/15_link_latency_ccdf.png", key="link_lat_ms",
                          title="Device-to-receiver link latency (excludes host stack/app)")
    if sdf is not None:
        plots.sensitivity(sdf, "cpu_load", [("interval_err_rms_us", "RMS interval error (µs)"),
                                            ("acq_sample_loss_rate", "acquisition loss rate"),
                                            ("lat_p99_ms", "p99 latency (ms)")],
                          f"{P}/16_sens_cpu_load.png", "CPU background load")
        plots.sensitivity(sdf, "buffer_size", [("rx_missing_rate", "receiver missing rate"),
                                               ("device_buffer_overflow_pkts", "device overflow (pkts)"),
                                               ("receiver_buffer_overflow_pkts", "receiver overflow (pkts)")],
                          f"{P}/17_sens_buffer.png", "device TX buffer (packets)", logx=True)
        plots.sensitivity(sdf, "stall_duration", [("rx_missing_rate", "receiver missing rate"),
                                                  ("lat_p99_ms", "p99 latency (ms)"),
                                                  ("acq_sample_loss_rate", "acquisition loss rate")],
                          f"{P}/18_sens_stall.png", "median stall duration (s)", logx=True)
        plots.sensitivity(sdf, "pio_ber", [("rx_missing_rate", "receiver missing rate"),
                                           ("crc_detected", "CRC-rejected packets")],
                          f"{P}/19_sens_ber.png", "PIO link bit-error rate", logx=True)

    # ---------------- manifest
    import matplotlib
    import scipy
    manifest = {
        "execution_date_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "runtime_s": round(time.time() - t0, 1),
        "python": sys.version.split()[0], "platform": platform.platform(),
        "numpy": np.__version__, "scipy": scipy.__version__, "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
        "base_seed": base["seed"],
        "seed_scheme": "np.random.SeedSequence([base_seed, trial]).spawn(5) -> streams "
                       "[ecg, measurement/ADC noise, CPU workload, ADC/DMA hardware, transport]; "
                       "identical across conditions and architectures (common random numbers)",
        "trial_indices": list(range(n_trials)),
        "n_trials_per_condition": n_trials, "duration_s": base["duration_s"], "fs": base["fs"],
        "conditions": base["conditions"], "quick_mode": args.quick,
        "data_provenance": "ALL results are simulation outputs; no hardware measurements were made.",
    }
    json.dump(manifest, open(os.path.join(out, "manifest.json"), "w"), indent=2)
    print(f"done in {time.time() - t0:.0f} s -> {out}")


if __name__ == "__main__":
    main()
