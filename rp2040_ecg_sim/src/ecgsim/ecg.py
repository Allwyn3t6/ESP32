"""Synthetic ECG generator and AD8232 / RP2040-ADC front-end model.

The ECG is a sum of Gaussian waves (P, Q, R, S, T) per beat, in the spirit of
McSharry et al. (2003), but evaluated directly in the time domain so that the
analog waveform can be sampled at *arbitrary* instants. That property is what
lets the simulator apply each architecture's real sampling instants to the same
continuous-time signal.
"""
from __future__ import annotations

import numpy as np

# (offset from R peak [s] at RR=1 s, amplitude [mV], width sigma [s]).
# P/T offsets are scaled by sqrt(RR) (Bazett-like); Q/R/S are not.
WAVES = {
    "P": (-0.200, 0.15, 0.025, True),
    "Q": (-0.035, -0.12, 0.010, False),
    "R": (0.000, 1.20, 0.010, False),
    "S": (0.035, -0.25, 0.010, False),
    "T": (0.280, 0.30, 0.050, True),
}


class SyntheticECG:
    """Continuous-time synthetic ECG (mV) with known R-peak times."""

    def __init__(self, cfg: dict, duration_s: float, rng: np.random.Generator):
        e = cfg["ecg"]
        self.cfg = e
        rr_mean = 60.0 / e["heart_rate_bpm"]
        n_beats = int(duration_s / rr_mean * 1.2) + 10
        # AR(1) RR variability, clipped to physiological range
        z = np.zeros(n_beats)
        eps = rng.normal(0, 1, n_beats)
        for i in range(1, n_beats):
            z[i] = 0.7 * z[i - 1] + np.sqrt(1 - 0.49) * eps[i]
        rr = np.clip(rr_mean * (1 + e["rr_sd_frac"] * z), 0.3, 2.0)
        beats = 0.35 + np.concatenate([[0.0], np.cumsum(rr[:-1])])
        self.rr = rr
        self.beats = beats
        self.r_times = beats[beats < duration_s - 0.3]
        self.duration_s = duration_s
        self.phi1, self.phi2, self.phi_m = rng.uniform(0, 2 * np.pi, 3)

    def clean(self, t: np.ndarray) -> np.ndarray:
        """Noise-free ECG + baseline wander in mV at times t."""
        t = np.asarray(t, dtype=float)
        idx = np.searchsorted(self.beats, t)
        out = np.zeros_like(t)
        nb = len(self.beats)
        for d in (-2, -1, 0, 1):
            jr = idx + d
            valid = (jr >= 0) & (jr < nb)
            j = np.clip(jr, 0, nb - 1)
            tb = np.where(valid, self.beats[j], np.inf)
            rr = self.rr[j]
            for off, amp, sig, scale in WAVES.values():
                o = off * np.sqrt(rr) if scale else off
                out += amp * np.exp(-0.5 * ((t - tb - o) / sig) ** 2)
        e = self.cfg
        out += e["wander_mv"] * np.sin(2 * np.pi * e["wander_hz"] * t + self.phi1)
        out += e["wander2_mv"] * np.sin(2 * np.pi * e["wander2_hz"] * t + self.phi2)
        return out

    def mains(self, t: np.ndarray) -> np.ndarray:
        e = self.cfg
        return e["mains_mv"] * np.sin(2 * np.pi * e["mains_hz"] * t + self.phi_m)


def afe_voltage(cfg: dict, ecg_mv: np.ndarray) -> np.ndarray:
    """AD8232 output model: V = offset + gain*ECG, clipped to the supply rails."""
    a = cfg["afe"]
    return np.clip(a["offset_v"] + a["gain_v_per_mv"] * ecg_mv, a["vmin"], a["vmax"])


def quantize(cfg: dict, v: np.ndarray, adc_noise_lsb: np.ndarray | None = None) -> np.ndarray:
    """Ideal 12-bit ADC (code = floor(V/Vref * 2^bits)), optional input-referred
    Gaussian noise in LSB, saturated to [0, 2^bits-1]. Returns uint16 codes."""
    a = cfg["adc"]
    full = 2 ** a["bits"]
    x = v / a["vref"] * full
    if adc_noise_lsb is not None:
        x = x + adc_noise_lsb
    return np.clip(np.floor(x), 0, full - 1).astype(np.uint16)


def code_to_volts(cfg: dict, code: np.ndarray) -> np.ndarray:
    """Mid-tread reconstruction of an ADC code to volts."""
    a = cfg["adc"]
    return (np.asarray(code, dtype=float) + 0.5) * a["vref"] / 2 ** a["bits"]
