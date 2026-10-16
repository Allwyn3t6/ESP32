import copy
import json
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from ecgsim import packets as P  # noqa: E402
from ecgsim.analysis import holm, paired_tests  # noqa: E402
from ecgsim.architectures import SharedAcquisitionDraws, deliver_bytes, run_architecture  # noqa: E402
from ecgsim.ecg import SyntheticECG, afe_voltage, code_to_volts, quantize  # noqa: E402
from ecgsim.experiment import apply_overrides, run_trial  # noqa: E402
from ecgsim.pipeline import Pipeline, SharedTransportDraws  # noqa: E402
from ecgsim.receiver import reconstruct  # noqa: E402
from ecgsim.rpeak import detect_r_peaks, score  # noqa: E402
from ecgsim.timing import Intervals  # noqa: E402


@pytest.fixture
def cfg():
    c = json.load(open(os.path.join(ROOT, "config.json")))
    c["duration_s"] = 30
    return c


def _ideal(cfg):
    return apply_overrides(cfg, cfg["conditions"]["baseline_no_faults"])


# ---------------- sample counts / quantisation / timestamps
def test_expected_number_of_samples(cfg):
    rows, _ = run_trial(_ideal(cfg), 0)
    n = cfg["fs"] * cfg["duration_s"]
    for r in rows:
        assert r["packets_sent"] + r["device_buffer_overflow_pkts"] == n // cfg["samples_per_packet"]
    # B and C: no faults -> every sample present
    assert rows[1]["rx_missing_rate"] == 0 and rows[2]["rx_missing_rate"] == 0


def test_adc_quantization_limits(cfg):
    v = np.array([-1.0, 0.0, 3.3 / 4096 * 0.999, 3.3 / 4096, 1.65, 3.2999, 3.3, 5.0])
    c = quantize(cfg, v)
    assert c.dtype == np.uint16
    assert c.min() >= 0 and c.max() <= 4095
    assert list(c) == [0, 0, 0, 1, 2048, 4095, 4095, 4095]
    # reconstruction error bounded by 1/2 LSB inside range
    x = np.linspace(0.001, 3.29, 10000)
    assert np.max(np.abs(code_to_volts(cfg, quantize(cfg, x)) - x)) <= 3.3 / 4096 / 2 + 1e-12
    # AFE clipping keeps the ADC input safe
    assert afe_voltage(cfg, np.array([-100.0, 100.0])).tolist() == [0.0, 3.3]


def test_timestamp_and_index_progression(cfg):
    c = _ideal(cfg)
    n = c["fs"] * c["duration_s"]
    ss = np.random.SeedSequence([1, 2]).spawn(5)
    rg = [np.random.default_rng(s) for s in ss]
    ecg = SyntheticECG(c, c["duration_s"], rg[0])
    acq = SharedAcquisitionDraws(c, rg[1], rg[2], rg[3], n, c["duration_s"])
    td = SharedTransportDraws(c, rg[4], n // c["samples_per_packet"], c["duration_s"])
    for arch in "ABC":
        res = run_architecture(arch, c, ecg, acq, td, n)
        t = res["t_acq"][~np.isnan(res["t_acq"])]
        assert np.all(np.diff(t) > 0), arch                      # strictly increasing instants
        ts = 1 / c["fs"]
        idx = np.nonzero(~np.isnan(res["t_acq"]))[0]
        assert np.all(res["t_acq"][idx] >= idx * ts - 1e-8)      # never sampled before its deadline (jitter aside)
        items, _ = deliver_bytes(res, c, td)
        seqs = [P.parse_packet(it[1])[0] for it in items]
        assert seqs == sorted(seqs) and np.all(np.diff(seqs) == c["samples_per_packet"])
        rx_t = [it[0] for it in items]
        assert np.all(np.diff(rx_t) >= 0)                        # receiver timestamps monotone
        assert np.all(np.array(rx_t) > res["ready"][: len(rx_t)])  # delivery after packet ready


# ---------------- packet / receiver
def _stream(n=500, N=50):
    codes = (np.arange(n) * 7 % 4096).astype(np.uint16)
    pk = [P.build_packet(s, codes[s:s + N]) for s in range(0, n, N)]
    return codes, pk


def test_packet_roundtrip_and_reconstruction():
    codes, pk = _stream()
    rx = reconstruct([(i, b) for i, b in enumerate(pk)], len(codes))
    assert np.array_equal(rx["codes"].astype(int), codes.astype(int))
    assert rx["duplicates"] == rx["out_of_order"] == rx["crc_fail"] == 0
    assert len(pk[0]) == P.packet_len(50) == 110


def test_detect_dropped_samples():
    codes, pk = _stream()
    del pk[3]
    rx = reconstruct([(i, b) for i, b in enumerate(pk)], len(codes))
    miss = np.isnan(rx["codes"])
    assert miss.sum() == 50 and miss[150:200].all()


def test_detect_acquisition_missing_marker():
    codes = np.arange(50, dtype=np.uint16)
    codes[[3, 10]] = P.MISSING
    rx = reconstruct([(0, P.build_packet(0, codes))], 50)
    assert rx["acq_marker"].sum() == 2 and np.isnan(rx["codes"]).sum() == 2 and rx["got"].all()


def test_detect_duplicates_and_reordering():
    codes, pk = _stream()
    items = [(0, pk[0]), (1, pk[2]), (2, pk[1]), (3, pk[1])] + [(i + 4, b) for i, b in enumerate(pk[3:])]
    rx = reconstruct(items, len(codes))
    assert rx["duplicates"] == 50
    assert rx["out_of_order"] == 50
    assert np.array_equal(rx["codes"].astype(int), codes.astype(int))


def test_crc_detects_corruption():
    _, pk = _stream()
    rng = np.random.default_rng(0)
    for _ in range(500):
        k = int(rng.integers(1, 4))
        bad = P.flip_bits(pk[0], rng.choice(len(pk[0]) * 8, k, replace=False))
        assert P.parse_packet(bad) is None   # CRC-16 detects all 1-3 bit errors at this length


# ---------------- buffers
def test_device_buffer_overflow_detection(cfg):
    c = apply_overrides(cfg, {"buffers.device_tx_packets": 2, "usb.stall_rate_hz": 0, "host.stall_rate_hz": 0})
    td = SharedTransportDraws(c, np.random.default_rng(0), 100, 10)
    td.usb_stalls = Intervals([0.05], [5.0])           # host stops polling for 5 s
    pipe = Pipeline(c, "usb", td, 110)
    for k in range(20):
        pipe.submit(k, 0.1 * k)
    # packet 0 (t=0) leaves before the stall at 50 ms; packets 1-2 fill the 2-slot buffer; 3..19 overflow
    assert pipe.dev_overflow == 17


def test_receiver_overflow_without_flow_control(cfg):
    c = apply_overrides(cfg, {"buffers.receiver_packets": 3, "host.stall_rate_hz": 0})
    td = SharedTransportDraws(c, np.random.default_rng(0), 100, 20)
    td.app_stalls = Intervals([0.0], [5.0])            # receiver app frozen for 5 s
    pio = Pipeline(c, "pio", td, 110)
    usb = Pipeline(c, "usb", td, 110)
    for k in range(40):
        pio.submit(k, 0.1 * k)
        usb.submit(k, 0.1 * k)
    assert pio.rx_overflow > 0                          # PIO: drops at receiver
    assert usb.rx_overflow == 0 and usb.dev_overflow > 0  # USB: back-pressure moves loss to device


def test_adc_fifo_overflow_detection(cfg):
    c = _ideal(cfg)
    c["adc"]["dma_latency_mean_s"] = 0.02             # pathological DMA starvation (10 sample periods)
    rows, _ = run_trial(c, 0)
    assert rows[1]["adc_fifo_overflow"] > 0 and rows[1]["acq_sample_loss_rate"] > 0


# ---------------- end-to-end consistency with configuration
def test_metrics_consistent_with_config(cfg):
    b = _ideal(cfg)
    rows, _ = run_trial(b, 0)
    A, B, C = rows
    assert B["interval_err_rms_us"] < 0.01 and C["interval_err_rms_us"] < 0.01   # 1 ns jitter configured
    assert A["interval_err_rms_us"] > B["interval_err_rms_us"]
    assert B["rmse_vs_afe_mv"] == pytest.approx(C["rmse_vs_afe_mv"])             # same acquisition path
    lsb_mv = 3300 / 4096
    expected = lsb_mv * np.sqrt(1 / 12 + b["adc"]["noise_lsb"] ** 2)
    assert B["rmse_vs_afe_mv"] == pytest.approx(expected, rel=0.1)
    for r in rows:
        assert r["r_sensitivity"] > 0.99 and r["r_ppv"] > 0.99
    # injected loss rate is recovered by sequence-gap detection
    lossy = apply_overrides(b, {"faults.loss_prob": 0.05})
    lossy["duration_s"] = 120
    rows, _ = run_trial(lossy, 3)
    for r in rows:
        assert 0.02 < r["packet_loss_rate"] < 0.09


def test_rpeak_scoring():
    true = np.array([1.0, 2.0, 3.0])
    tp, fp, fn, e = score(true, np.array([1.01, 2.2, 3.0, 4.0]))
    assert (tp, fp, fn) == (2, 2, 1) and np.allclose(sorted(np.abs(e)), [0.0, 0.01])


def test_rpeak_detector_on_clean_ecg(cfg):
    ecg = SyntheticECG(cfg, 30, np.random.default_rng(5))
    t = np.arange(30 * 500) / 500
    det = detect_r_peaks(afe_voltage(cfg, ecg.clean(t)), 500)
    det = det[(det > 0.2) & (det < 29.7)]
    tp, fp, fn, e = score(ecg.r_times, det)
    assert fn == 0 and fp == 0 and np.sqrt(np.mean(e ** 2)) < 2e-3


# ---------------- reproducibility
def test_reproducible_with_fixed_seed(cfg):
    c = copy.deepcopy(cfg)
    r1, _ = run_trial(c, 7)
    r2, _ = run_trial(c, 7)
    r3, _ = run_trial(c, 8)
    assert r1 == r2 or all(
        all((a[k] == b[k]) or (isinstance(a[k], float) and np.isnan(a[k]) and np.isnan(b[k])) for k in a)
        for a, b in zip(r1, r2))
    assert r1[0]["interval_err_rms_us"] != r3[0]["interval_err_rms_us"]


def test_holm_and_paired_tests():
    p = holm(np.array([0.01, 0.04, 0.03]))
    assert np.allclose(p, [0.03, 0.06, 0.06])
    import pandas as pd
    rows = []
    for t in range(10):
        for a, v in (("A", 10 + t), ("B", 1.0), ("C", 1.0)):
            rows.append({"condition": "x", "trial": t, "arch": a, **{m: v for m in
                         ["interval_err_rms_us", "acq_sample_loss_rate", "rx_missing_rate", "packet_loss_rate",
                          "device_buffer_overflow_pkts", "receiver_buffer_overflow_pkts", "lat_mean_ms",
                          "lat_p99_ms", "lat_jitter_std_ms", "link_lat_jitter_std_ms", "rmse_vs_afe_mv",
                          "rmse_filled_mv", "r_sensitivity", "r_timing_err_rms_ms"]}})
    out = paired_tests(pd.DataFrame(rows))
    bc = out[out.comparison == "B vs C"]
    assert (bc.better == "no significant difference").all()
