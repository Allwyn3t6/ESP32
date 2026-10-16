# Hardware validation roadmap (Raspberry Pi Pico + AD8232)

**Status: future work.** Nothing in this document has been built or measured. The simulation
(`REPORT.md`) produces *predictions*; this plan describes how to obtain *hardware evidence* that
confirms or refutes them. Keep the two kinds of evidence separate in any publication.

## 0. Safety

ECG electrodes are connected to a person. During development, power the Pico and AD8232 from a
battery or an isolated USB supply / USB isolator (e.g. ADuM4160-based) and never connect a
mains-powered PC to a subject without medical-grade isolation. Use an ECG simulator/patient
simulator or a function generator with an attenuator whenever possible.

## 1. Bench signal sources (in order of preference)

1. Arbitrary waveform generator replaying the *same* synthetic ECG used in the simulation
   (export `SyntheticECG.clean()` at ≥ 50 kS/s), attenuated to mV level into the AD8232 inputs,
   so received data can be compared to a known reference exactly as in the simulation.
2. Commercial ECG patient simulator (known heart rate).
3. Human subject (only for qualitative checks; no ground truth).

## 2. Firmware (Pico SDK, C)

Start from the official example `adc/dma_capture/dma_capture.c` in
https://github.com/raspberrypi/pico-examples (ADC free-running, `adc_fifo_setup(true, true, 1, false, false)`,
`adc_set_clkdiv()`, DMA paced by `DREQ_ADC`). Build three firmware images sharing one acquisition module.

### Common
* `adc_init(); adc_gpio_init(26); adc_select_input(0);` (AD8232 OUTPUT → GPIO26/ADC0, AD8232 powered from 3V3).
* Sample word: 12-bit code in a `uint16_t`; enable the FIFO error bit (`adc_fifo_setup(..., err_in_fifo=true, ...)`)
  so conversion errors are visible; `0xFFFF` marks a missed slot (as in the simulation).
* Packet format identical to `src/ecgsim/packets.py` (magic, seq32, count, flags, N×u16, CRC-16/CCITT-FALSE).
  For C, compute the CRC with the DMA sniffer (`dma_sniffer_enable(chan, DMA_SNIFF_CTRL_CALC_VALUE_CRC16, true)`)
  and verify it matches the host implementation.
* Instrumentation GPIOs (logic analyser): `GP2` toggled at each acquisition (A: right after `adc_read()`;
  B/C: driven from a PWM slice or ADC-synchronous PIO SM, or toggle in the DMA IRQ and log per-block timing),
  `GP3` high while a packet is handed to the transport, `GP4` asserted on any overflow event.

### A — CPU polling
```c
uint64_t next = time_us_64();
for (;;) {
    if (time_us_64() >= next) { gpio_xor_mask(1u<<2); slot[i++] = adc_read(); next += 2000; /* skip-ahead if late */ }
    if (i == N) { build_packet(); fwrite(pkt, 1, len, stdout); i = 0; }   /* stdio_usb, default timeout */
    background_tasks(); tud_task();
}
```
Record missed deadlines (when `time_us_64() - next >= 2000`).

### B — ADC FIFO + DMA ping-pong + USB CDC
* `adc_set_clkdiv(95999)` → 48 MHz / 96 000 = 500 S/s exactly (integer divider: no fractional dither).
* Two DMA channels chained to each other (`channel_config_set_chain_to`), each writing N samples into its own buffer,
  `DREQ_ADC`, 16-bit transfers, `adc_fifo_setup(..., dreq_en=true, dreq_thresh=1, ...)`.
* DMA IRQ: mark buffer ready, re-arm its write address; main loop: `tud_cdc_write(pkt) ; tud_cdc_write_flush()`
  and `tud_task()`. Count "buffer not consumed before DMA wrapped" as DMA-overflow; read `ADC_FCS.OVER/UNDER`.

### C — B + PIO transmitter
* Same acquisition as B.
* PIO program: synchronous serial (data + clock) or UART-style 8N1 TX at 2 Mbit/s; optional `CTS` input pin for
  the flow-control variant suggested by the simulation (`wait 1 pin` before each byte).
* DMA channel 3: memory → `pio->txf[sm]`, `DREQ_PIO0_TX0`, triggered (chained) when a ping-pong half completes,
  so packet emission needs no CPU except writing the 8-byte header.
* Receiver: FT232H in synchronous FIFO mode, a second Pico running a PIO receiver with USB to the PC, or a
  logic analyser with protocol decode. The receiver must hardware-timestamp packet arrival.

## 3. Host receiver (Python)

* USB CDC: `pyserial` read loop in a dedicated thread; timestamp with `time.perf_counter_ns()` on read return.
* PIO link: receiver-side timestamps + PC arrival timestamps.
* Reuse `src/ecgsim/receiver.py::reconstruct` unchanged so the metrics are computed identically to the simulation.

## 4. Instrumentation and what to measure

| Simulation quantity | Hardware measurement |
|---|---|
| sampling-interval error | logic analyser (≥ 100 MS/s) on `GP2`; histogram of edge-to-edge intervals over ≥ 10 min |
| ADC conversion instant vs. deadline (A) | `GP2` vs a PWM reference running from the same clock |
| acquisition sample loss | missed-deadline counter (A), `ADC_FCS.OVER` (B/C) |
| device buffer overflow | firmware counters + `GP4` |
| link latency | `GP3` (packet handed to transport) → receiver hardware timestamp (FT232H / second Pico / USB protocol analyser for CDC) |
| end-to-end latency | `GP3` edge captured by a host-synchronised input, or loop-back: host timestamps of a GPIO-triggered marker packet |
| packet loss / dup / reorder | sequence numbers (receiver.py) |
| corruption | CRC failures at receiver; inject known bit errors on the PIO line with a series switch/noise source to verify detection |
| CPU utilisation | toggle a GPIO while in acquisition/transport code; duty cycle on a logic analyser; or `time_us_64()` accumulators |
| waveform RMSE, R-peaks | AWG replaying the known synthetic ECG; align by cross-correlation; reuse `rpeak.py` |

Stress conditions to replicate: synthetic CPU load (busy loops of lognormal duration in the main loop, matched to
`config.json` `cpu.*`), host load (`stress-ng` on the PC), forced read stalls in the host script, reduced TinyUSB
buffer sizes (`CFG_TUD_CDC_TX_BUFSIZE`), long/noisy cables for the PIO link.

## 5. Decision criteria the hardware must settle

The simulation could not settle these; they require measurements:
1. Real distribution of PC-side USB CDC delivery latency and stall episodes on the target host OS (the dominant term in
   every latency result).
2. Real bit-error rate of the PIO link for the intended cable, rate and environment (decides whether C needs FEC/ARQ).
3. Actual CPU workload of the intended application (filtering, display, BLE...) → A's missed-sample rate.
4. RP2040 ADC effective noise and DNL with the AD8232 connected (sets the waveform-error floor for all architectures).
5. Whether a flow-control line on the PIO link (or a sufficiently large receiver buffer) removes C's receiver-overflow
   failure mode while keeping its deterministic link timing.

Pre-register the comparison: same AWG input, same firmware acquisition module for B and C, ≥ 30 ten-minute runs per
condition, paired analysis with the scripts in this repository.
