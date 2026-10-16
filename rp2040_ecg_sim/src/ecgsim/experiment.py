"""One Monte Carlo trial = one ECG record pushed through architectures A, B and C."""
from __future__ import annotations

import copy

import numpy as np

from .architectures import SharedAcquisitionDraws, deliver_bytes, run_architecture
from .ecg import SyntheticECG, afe_voltage, code_to_volts
from .pipeline import SharedTransportDraws
from .receiver import reconstruct
from .rpeak import detect_r_peaks, score

ARCHS = ("A", "B", "C")
ARCH_NAMES = {"A": "A: CPU polling + USB CDC", "B": "B: ADC FIFO+DMA + USB CDC", "C": "C: ADC FIFO+DMA + PIO/DMA link"}


def apply_overrides(cfg: dict, overrides: dict) -> dict:
    cfg = copy.deepcopy(cfg)
    for key, val in overrides.items():
        if key.startswith("_"):
            continue
        sec, name = key.split(".")
        if sec not in cfg or name not in cfg[sec]:
            raise KeyError(f"unknown config key {key}")
        cfg[sec][name] = val
    return cfg


def trial_seed_sequence(base_seed: int, trial: int) -> np.random.SeedSequence:
    """Seeds depend only on (base_seed, trial): conditions/architectures share them (CRN)."""
    return np.random.SeedSequence([base_seed, trial])


def _pct(x, q):
    return float(np.percentile(x, q)) if len(x) else np.nan


def cpu_utilisation_estimate(cfg, arch):
    """Derived (not measured) CPU-time fraction spent on acquisition+transmission
    from the assumed per-operation costs in config['cpu']."""
    c = cfg["cpu"]
    pkt_period = cfg["samples_per_packet"] / cfg["fs"]
    if arch == "A":
        return 1.0   # busy-wait polling loop occupies the core continuously
    if arch == "B":
        return (c["irq_cost_s"] + c["usb_task_cost_s"]) / pkt_period
    return c["irq_cost_s"] / pkt_period


def run_trial(cfg: dict, trial: int, keep_raw: bool = False):
    fs, T = cfg["fs"], cfg["duration_s"]
    n = int(round(fs * T))
    N = cfg["samples_per_packet"]
    n_pk = int(np.ceil(n / N))
    Ts = 1.0 / fs
    ss = trial_seed_sequence(cfg["seed"], trial)
    r_ecg, r_noise, r_cpu, r_hw, r_tr = [np.random.default_rng(s) for s in ss.spawn(5)]
    ecg = SyntheticECG(cfg, T, r_ecg)
    acq = SharedAcquisitionDraws(cfg, r_noise, r_cpu, r_hw, n, T)
    tdraws = SharedTransportDraws(cfg, r_tr, n_pk, T)
    t_nom = np.arange(n) * Ts
    ref_v = afe_voltage(cfg, ecg.clean(t_nom))          # noise-free reference at ideal instants
    ref_ac = ref_v - ref_v.mean()
    # AFE output incl. measurement noise at ideal instants: what a perfect sampler would see
    ref_in = afe_voltage(cfg, ecg.clean(t_nom) + ecg.mains(t_nom) + acq.white)
    rows, raw = [], {}
    for arch in ARCHS:
        res = run_architecture(arch, cfg, ecg, acq, tdraws, n)
        items, n_corrupt = deliver_bytes(res, cfg, tdraws)
        rx = reconstruct(items, n)
        m = {"trial": trial, "arch": arch, "seed_entropy": cfg["seed"], "seed_spawn_key": trial}
        # --- acquisition timing (ground-truth instants, independent of host timing)
        ta = res["t_acq"]
        ok = ~np.isnan(ta)
        idx = np.nonzero(ok)[0]
        dI = np.diff(ta[idx]) - np.diff(idx) * Ts
        m["interval_err_rms_us"] = float(np.sqrt(np.mean(dI ** 2)) * 1e6)
        m["interval_err_std_us"] = float(np.std(dI) * 1e6)
        m["interval_err_maxabs_us"] = float(np.max(np.abs(dI)) * 1e6)
        m["abs_time_err_rms_us"] = float(np.sqrt(np.mean((ta[idx] - t_nom[idx]) ** 2)) * 1e6)
        m["acq_sample_loss_rate"] = float(res["acq_missed"].mean())
        m["adc_fifo_overflow"] = int(res["fifo_overflow"].sum())
        m["device_buffer_overflow_pkts"] = res["pipe"].dev_overflow
        m["receiver_buffer_overflow_pkts"] = res["pipe"].rx_overflow
        # --- receiver integrity
        have = ~np.isnan(rx["codes"])
        m["rx_missing_rate"] = float(1 - have.mean())
        m["rx_transport_missing_rate"] = float((~rx["got"]).mean())
        m["duplicates"] = rx["duplicates"]
        m["out_of_order"] = rx["out_of_order"]
        slices = res["slices"]
        pk_absent = sum(1 for a0, a1 in slices if not rx["got"][a0:a1].any())
        m["packets_sent"] = len(slices) - res["pipe"].dev_overflow
        m["packets_lost_detected"] = pk_absent
        m["packet_loss_rate"] = pk_absent / len(slices)
        m["corrupted_packets"] = n_corrupt
        m["crc_detected"] = rx["crc_fail"]
        m["corruption_detection_rate"] = rx["crc_fail"] / n_corrupt if n_corrupt else np.nan
        # --- latency (first delivered copy of each packet)
        first = {}
        link = {}
        for it in items:
            k = it[2]
            if k not in first:
                first[k] = it[0]
                link[k] = it[3]
        ks = np.array(sorted(first))
        if ks.size:
            lat = np.array([first[k] for k in ks]) - res["ready"][ks]
            llat = np.array([link[k] for k in ks]) - res["ready"][ks]
        else:
            lat = llat = np.array([])
        lat_ms, llat_ms = lat * 1e3, llat * 1e3
        m["lat_mean_ms"] = float(np.mean(lat_ms)) if lat.size else np.nan
        m["lat_median_ms"] = _pct(lat_ms, 50)
        m["lat_p95_ms"] = _pct(lat_ms, 95)
        m["lat_p99_ms"] = _pct(lat_ms, 99)
        m["lat_max_ms"] = float(np.max(lat_ms)) if lat.size else np.nan
        m["lat_jitter_std_ms"] = float(np.std(lat_ms)) if lat.size else np.nan
        m["link_lat_mean_ms"] = float(np.mean(llat_ms)) if lat.size else np.nan
        m["link_lat_p99_ms"] = _pct(llat_ms, 99)
        m["link_lat_jitter_std_ms"] = float(np.std(llat_ms)) if lat.size else np.nan
        # --- waveform
        v = code_to_volts(cfg, rx["codes"])
        err = v[have] - ref_v[have]
        m["rmse_mv"] = float(np.sqrt(np.mean(err ** 2)) * 1e3)
        m["snr_db"] = float(10 * np.log10(np.mean(ref_ac[have] ** 2) / np.mean(err ** 2)))
        m["rmse_vs_afe_mv"] = float(np.sqrt(np.mean((v[have] - ref_in[have]) ** 2)) * 1e3)
        xi = np.arange(n)
        vf = np.interp(xi, xi[have], v[have])
        m["rmse_filled_mv"] = float(np.sqrt(np.mean((vf - ref_v) ** 2)) * 1e3)
        # --- R peaks
        det = detect_r_peaks(vf, fs)
        det = det[(det > 0.2) & (det < T - 0.3)]          # same scoring window as the truth set
        tp, fp, fn, e = score(ecg.r_times, det)
        m["r_true"] = len(ecg.r_times)
        m["r_tp"], m["r_fp"], m["r_fn"] = tp, fp, fn
        m["r_sensitivity"] = tp / (tp + fn) if tp + fn else np.nan
        m["r_ppv"] = tp / (tp + fp) if tp + fp else np.nan
        m["r_timing_err_rms_ms"] = float(np.sqrt(np.mean(e ** 2)) * 1e3) if e.size else np.nan
        m["r_timing_err_maxabs_ms"] = float(np.max(np.abs(e)) * 1e3) if e.size else np.nan
        m["cpu_util_est"] = cpu_utilisation_estimate(cfg, arch)
        rows.append(m)
        if keep_raw:
            raw[arch] = dict(lat_ms=lat_ms, link_lat_ms=llat_ms, interval_err_us=dI * 1e6,
                             v=v[: 10 * fs].copy(), ref=ref_v[: 10 * fs].copy(), r_err_ms=e * 1e3)
    return rows, raw
