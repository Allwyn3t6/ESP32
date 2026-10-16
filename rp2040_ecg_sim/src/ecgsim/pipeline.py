"""Packet-level discrete-event model of: device TX buffer -> link -> receiver buffer -> receiver app.

Links
-----
* "usb": USB CDC (full-speed bulk IN). Host-polled; flow controlled (device cannot
  send when the host-side buffer is full, it NAKs instead). Line errors are handled
  by USB's link-layer CRC + retry, so they appear as (negligible) delay, not corruption.
  Subject to host-controller/OS polling stalls.
* "pio": unidirectional PIO-driven serial link at a fixed bit rate, fed by DMA.
  No flow control (no handshake line), so a full receiver buffer drops packets.
  Line bit errors (BER) corrupt packets; the receiver's CRC detects them.

The packets are processed in order, one at a time; all state needed for packet k
depends only on packets < k, so a sequential loop is an exact event simulation of
this queueing network.
"""
from __future__ import annotations

from collections import deque

import numpy as np

from .timing import Intervals, lognormal_median


class SharedTransportDraws:
    """Per-packet random draws shared by all architectures (common random numbers)."""

    def __init__(self, cfg, rng, n_packets, T):
        h, u = cfg["host"], cfg["usb"]
        self.host_lat = lognormal_median(rng, h["latency_median_s"], h["latency_sigma"], n_packets)
        self.proc = lognormal_median(rng, h["proc_median_s"], h["proc_sigma"], n_packets)
        self.u_frame = rng.uniform(0, 1, n_packets)
        self.rx_jit = rng.normal(0, 1, n_packets)
        self.usb_stalls = Intervals.poisson(
            rng, u["stall_rate_hz"], T + 5,
            lambda n: lognormal_median(rng, u["stall_median_s"], u["stall_sigma"], n), queue=False)
        self.app_stalls = Intervals.poisson(
            rng, h["stall_rate_hz"], T + 5,
            lambda n: lognormal_median(rng, h["stall_median_s"], h["stall_sigma"], n), queue=False)
        f = cfg["faults"]
        self.lost = rng.uniform(size=n_packets) < f["loss_prob"]
        self.corrupt = rng.uniform(size=n_packets) < f["corrupt_prob"]
        self.dup = rng.uniform(size=n_packets) < f["dup_prob"]
        self.reorder = rng.uniform(size=n_packets) < f["reorder_prob"]
        self.n_inj_bits = rng.integers(1, 4, n_packets)
        self.bit_seed = rng.integers(0, 2 ** 31, n_packets)
        self.ber_u = rng.uniform(size=n_packets)


class Pipeline:
    def __init__(self, cfg, link: str, draws: SharedTransportDraws, pkt_bytes: int):
        assert link in ("usb", "pio")
        self.cfg, self.link, self.d = cfg, link, draws
        self.qdev = int(cfg["buffers"]["device_tx_packets"])
        self.qrx = int(cfg["buffers"]["receiver_packets"])
        self.pkt_bytes = pkt_bytes
        if link == "usb":
            u = cfg["usb"]
            self.tx_time = pkt_bytes * 8 * u["protocol_overhead"] / u["bitrate"]
        else:
            p = cfg["pio"]
            self.tx_time = pkt_bytes * p["bits_per_byte"] / p["bitrate"]
        self.dev_q = deque()          # departure (tx start) times of packets in device buffer
        self.last_tx_end = 0.0
        self.last_arrival = 0.0
        self.last_finish = 0.0
        self.rx_app_starts = deque()  # app start times of packets accepted into receiver buffer
        self.rx_all_starts = []       # for flow control indexing
        self.dev_overflow = 0
        self.rx_overflow = 0
        self.delivered = []           # (k, deliver_time, link_arrival_time, is_dup)

    def submit(self, k: int, h: float, blocking: bool = False, timeout: float = 0.0):
        """Offer packet k to the device TX buffer at time h.

        Returns (accepted: bool, cpu_free_time). For blocking (CPU-driven) firmware
        cpu_free_time is when the writer returns; otherwise it equals h."""
        q = self.dev_q
        while q and q[0] <= h:
            q.popleft()
        t_acc = h
        if len(q) >= self.qdev:
            if blocking and q[0] - h <= timeout:
                t_acc = q[0]
                q.popleft()
            else:
                self.dev_overflow += 1
                return False, (h + timeout if blocking else h)
        # ---- link transmission
        d = self.d
        if self.link == "usb":
            start = max(t_acc, self.last_tx_end)
            if len(self.rx_all_starts) >= self.qrx:       # host buffer full -> NAK (flow control)
                start = max(start, self.rx_all_starts[-self.qrx])
            start += d.u_frame[k] * self.cfg["usb"]["frame_s"]
            start = d.usb_stalls.push(start)
            end = start + self.tx_time
            arrival_link = end
        else:
            start = max(t_acc, self.last_tx_end)
            end = start + self.tx_time
            p = self.cfg["pio"]
            arrival_link = end + p["rx_iface_latency_s"] + abs(d.rx_jit[k]) * p["rx_iface_jitter_s"]
        self.last_tx_end = end
        q.append(start)
        # ---- receiver side
        if d.lost[k]:
            return True, t_acc
        a = max(arrival_link + d.host_lat[k], self.last_arrival)
        self.last_arrival = a
        rq = self.rx_app_starts
        while rq and rq[0] <= a:
            rq.popleft()
        if self.link == "pio" and len(rq) >= self.qrx:
            self.rx_overflow += 1                       # no flow control -> drop
            return True, t_acc
        copies = 2 if d.dup[k] else 1
        for c in range(copies):
            s_app = d.app_stalls.push(max(a, self.last_finish))
            fin = s_app + d.proc[k]
            self.last_finish = fin
            rq.append(s_app)
            self.rx_all_starts.append(s_app)
            self.delivered.append((k, fin, arrival_link, c == 1))
        return True, t_acc
