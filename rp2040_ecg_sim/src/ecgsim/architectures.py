"""Acquisition + packetisation models for architectures A, B and C.

A  CPU polling superloop: timer-polled adc_read(); background tasks run to completion
   between deadline checks; ISRs preempt; packets written with a blocking stdio-like
   USB CDC write that gives up after a timeout. Missed deadlines are skipped
   (slot transmitted as 0xFFFF).
B  ADC free-running (CLKDIV) -> 4-deep FIFO -> DMA ping-pong; half-buffer IRQ queues a
   packet that the USB stack (tud_task in the superloop) sends via USB CDC.
C  Same acquisition as B; completed half-buffer is chained by DMA to a PIO state machine
   that serialises it on a custom link. CPU only builds the 8-byte header in the IRQ.
"""
from __future__ import annotations

import numpy as np

from . import packets as P
from .ecg import afe_voltage, quantize
from .pipeline import Pipeline
from .timing import Intervals, lognormal_mean


class SharedAcquisitionDraws:
    """Random quantities shared across architectures (common random numbers)."""

    def __init__(self, cfg, rng_noise, rng_cpu, rng_hw, n, T):
        c, a = cfg["cpu"], cfg["adc"]
        self.white = rng_noise.normal(0, cfg["ecg"]["noise_mv"], n)
        self.adc_noise = rng_noise.normal(0, a["noise_lsb"], n) if a["noise_lsb"] > 0 else None
        # CPU background workload (non-preemptive tasks, Poisson arrivals, utilisation = load)
        if c["load"] > 0:
            lam = c["load"] / c["task_mean_s"]
            m = rng_cpu.poisson(lam * (T + 1))
            arr = np.sort(rng_cpu.uniform(0, T + 1, m))
            dur = lognormal_mean(rng_cpu, c["task_mean_s"], c["task_sigma"], m)
            self.tasks = Intervals.fifo_queue(arr, dur)      # NOT merged: loop checks between tasks
        else:
            self.tasks = Intervals()
        self.isr = Intervals.poisson(rng_cpu, c["isr_rate_hz"], T + 1,
                                     lambda k: rng_cpu.uniform(c["isr_min_s"], c["isr_max_s"], k), queue=True)
        self.poll = rng_cpu.uniform(0, c["poll_granularity_s"], n) + rng_cpu.uniform(0, c["timer_resolution_s"], n)
        self.hw_jitter = rng_hw.normal(0, a["hw_jitter_rms_s"], n)
        self.dma_lat = rng_hw.exponential(a["dma_latency_mean_s"], n)


def sample_codes(cfg, ecg, t_sample, idx, draws):
    """Analog ECG evaluated at the true sampling instants -> AFE -> ADC codes."""
    mv = ecg.clean(t_sample) + ecg.mains(t_sample) + draws.white[idx]
    v = afe_voltage(cfg, mv)
    noise = draws.adc_noise[idx] if draws.adc_noise is not None else None
    return quantize(cfg, v, noise)


def _packet_slices(n, N):
    return [(s, min(s + N, n)) for s in range(0, n, N)]


def run_architecture(arch, cfg, ecg, acq: SharedAcquisitionDraws, tdraws, n):
    fs, N = cfg["fs"], cfg["samples_per_packet"]
    Ts = 1.0 / fs
    c = cfg["cpu"]
    deadlines = np.arange(n) * Ts
    slices = _packet_slices(n, N)
    pipe = Pipeline(cfg, "pio" if arch == "C" else "usb", tdraws, P.packet_len(N))
    t_acq = np.full(n, np.nan)          # ground-truth acquisition instants
    acq_missed = np.zeros(n, bool)
    fifo_overflow = np.zeros(n, bool)
    ready = np.zeros(len(slices))
    pkts = []

    if arch == "A":
        free = acq.tasks.end_if_inside(deadlines)
        base_missed = np.zeros(n, bool)
        base_missed[:-1] = free[:-1] >= deadlines[1:]
        s_base = acq.isr.end_if_inside(free) + acq.poll
        nxt = np.append(deadlines[1:], np.inf)
        u = 0.0                                   # CPU blocked in USB write until u
        for k, (a0, a1) in enumerate(slices):
            dl, nx = deadlines[a0:a1], nxt[a0:a1]
            s = s_base[a0:a1].copy()
            miss = base_missed[a0:a1].copy()
            blk = dl < u
            if blk.any():
                miss |= blk & (nx < u)            # deadline and its successor both elapsed while blocked
                take = blk & ~(nx < u)            # latest elapsed deadline sampled when write returns
                s[take] = np.maximum(s[take], u + acq.poll[a0:a1][take])
                miss |= take & (s >= nx)
            t_acq[a0:a1] = np.where(miss, np.nan, s)
            acq_missed[a0:a1] = miss
            last = np.nanmax(t_acq[a0:a1]) if (~miss).any() else max(u, dl[-1])
            h = max(last, u) + c["format_cost_s"]
            ready[k] = h
            codes = np.full(a1 - a0, P.MISSING, np.uint16)
            ok = ~miss
            if ok.any():
                codes[ok] = sample_codes(cfg, ecg, s[ok], np.arange(a0, a1)[ok], acq)
            pkts.append(P.build_packet(a0, codes))
            _, u = pipe.submit(k, h, blocking=True, timeout=c["stdio_timeout_s"])
            u = max(u, h)
    else:
        a = cfg["adc"]
        t_hw = deadlines + acq.hw_jitter
        conv_done = t_hw + a["conv_time_s"]
        rd = np.maximum.accumulate(conv_done + acq.dma_lat)   # DMA reads FIFO in order
        depth = a["fifo_depth"]
        fifo_overflow[depth:] = rd[:-depth] > conv_done[depth:]
        t_acq[:] = np.where(fifo_overflow, np.nan, t_hw)
        acq_missed[:] = fifo_overflow
        codes_all = np.full(n, P.MISSING, np.uint16)
        ok = ~fifo_overflow
        codes_all[ok] = sample_codes(cfg, ecg, t_hw[ok], np.nonzero(ok)[0], acq)
        for k, (a0, a1) in enumerate(slices):
            r = rd[a1 - 1]                                     # DMA half-buffer complete
            if arch == "B":
                h = float(acq.tasks.end_if_inside(r + c["irq_cost_s"])) + c["usb_task_cost_s"]
            else:
                h = r + c["irq_cost_s"] + cfg["pio"]["dma_chain_latency_s"]
            ready[k] = r
            pkts.append(P.build_packet(a0, codes_all[a0:a1]))
            pipe.submit(k, h, blocking=False)
    return dict(pipe=pipe, packets=pkts, ready=ready, t_acq=t_acq, acq_missed=acq_missed,
                fifo_overflow=fifo_overflow, slices=slices, link="pio" if arch == "C" else "usb")


def deliver_bytes(res, cfg, tdraws):
    """Apply line/host faults to delivered packets; return list of (deliver_time, bytes, k, meta)."""
    pipe, pkts = res["pipe"], res["packets"]
    out = []
    corrupted_true = 0
    ber = cfg["pio"]["ber"] if res["link"] == "pio" else 0.0
    for (k, t_del, t_link, is_dup) in pipe.delivered:
        b = pkts[k]
        nbits = len(b) * 8
        flips = []
        rng = np.random.default_rng(int(tdraws.bit_seed[k]) + (7919 if is_dup else 0))
        if ber > 0:
            # P(at least one error) test with a shared uniform, then exact binomial count
            p_any = 1 - (1 - ber) ** nbits
            if tdraws.ber_u[k] < p_any:
                m = max(1, rng.binomial(nbits, ber))
                flips += list(rng.choice(nbits, m, replace=False))
        if tdraws.corrupt[k] and not is_dup:
            flips += list(rng.choice(nbits, int(tdraws.n_inj_bits[k]), replace=False))
        if flips:
            b = P.flip_bits(b, flips)
            corrupted_true += 1
        out.append([t_del, b, k, t_link, bool(flips), is_dup])
    # host-side reordering: swap delivery order with the following item
    i = 0
    while i < len(out) - 1:
        if tdraws.reorder[out[i][2]] and not out[i][5]:
            out[i][1:], out[i + 1][1:] = out[i + 1][1:], out[i][1:]
            i += 2
        else:
            i += 1
    return out, corrupted_true
