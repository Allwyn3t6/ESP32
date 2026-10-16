# ECG acquisition on the RP2040: CPU polling vs ADC FIFO + DMA vs PIO + DMA — a Monte Carlo simulation study

**Evidence type: simulation only.** Every number in this report comes from the simulator in this repository
(executed 2026-10-16, see `results/manifest.json`). No Raspberry Pi Pico or AD8232 was measured. Where a parameter
is a documented hardware specification it is marked **[spec]**; everything else is a **[assumption]** or a
**[derived]** quantity. Full per-condition tables with 95 % CIs: `results/summary_tables.md`.

---

## Abstract

We simulated three ways of getting an AD8232 ECG signal from an RP2040 to a receiver: (A) a CPU superloop that polls
a timer and calls `adc_read()`, then writes packets to USB CDC; (B) the ADC free-running from its clock divider into
its FIFO, with DMA ping-pong buffers and USB CDC; (C) the same ADC + DMA acquisition as B, with a DMA-fed PIO state
machine that drives a custom unidirectional serial link. The same synthetic ECG (P/QRS/T, HRV, baseline wander, noise),
500 S/s sampling, 12-bit quantisation over 3.3 V and binary packet format were used for all three. We ran 8 conditions
× 30 Monte Carlo trials of 10-minute records (240 trials, 720 architecture-records, 216 million simulated samples), plus four sensitivity sweeps.

Sample timing in A was the clear weak point. Its RMS sampling-interval error was 260 µs under normal load and 723 µs
under high load. B and C sample with the ADC's hardware pacing and showed 1.4 ns of error, which is just the jitter we
assumed. Under high CPU load A missed 22.3 % of its samples; B and C missed none. **B and C have exactly the same
acquisition metrics**, as they must. The PIO stage made a statistically significant difference only in the
*digital transmission stage*:

- It cut link-latency jitter from 0.39 ms to 0.3 µs.
- It kept end-to-end p99 latency flat under CPU load (4.0 ms vs 15.7 ms for B), because it does not depend on the
  USB stack being serviced in the main loop.
- It lost fewer samples than B when the USB host stalled in bursts (0.007 % vs 0.198 %).

It was not better in every case:

- With small buffers, C dropped data at the receiver (no flow control) instead of at the device, and lost about the
  same amount as B (0.232 % vs 0.222 %, not significant).
- With a raised line bit-error rate, C lost more than B (1.18 % vs 1.11 %, p_Holm = 3·10⁻⁴), because USB corrects
  line errors and the PIO link only detects them.

The waveform and R-peak metrics did not differ between B and C except where packets were lost. **Answer to the
research question:** under this model, PIO + DMA gives a measurable improvement in *transmission timing determinism*
and in latency robustness to CPU load. It does not give a general improvement in *data integrity*. Whether the timing
advantage matters, and whether the integrity trade-off favours C, depends on parameters that only hardware
measurements can establish (§11).

---

## 1. Research objectives

1. Quantify the sampling determinism, data continuity, latency and signal fidelity of A, B and C under identical ECG
   input and quantisation.
2. Separate **acquisition timing** (when the ADC actually samples) from **transport timing** (when data reaches the
   receiver).
3. Test whether a PIO + DMA transmission stage (C) improves on ADC FIFO + DMA + USB CDC (B). Because B and C share
   the acquisition path, any difference must come from the transmission stage.
4. Find out which conclusions hold up when the assumptions change (sensitivity analysis), and list the measurements
   needed where they do not.

## 2. Technical background

- **RP2040 ADC [spec]** (RP2040 datasheet §4.9):
  - SAR ADC, 12 bits, clocked at 48 MHz; 96 clock cycles per conversion, i.e. 2 µs, so at most 500 kS/s.
  - In free-running mode the conversion interval is set by `DIV` as (1 + DIV_INT + DIV_FRAC/256) ADC clocks.
    500 S/s needs DIV_INT = 95 999 and a zero fraction, so there is no fractional-divider dither.
  - Results go into a 4-entry FIFO that can raise a DMA request (DREQ_ADC). FIFO overflow is flagged in `FCS.OVER`.
  - The datasheet reports an effective resolution well below 12 bits and an ADC DNL erratum (RP2040-E11). We model
    this as additive noise (§7).
- **DMA [spec]**: 12 channels that can be chained, paced by DREQs. A *sniffer* can compute a CRC (including
  CRC-16-CCITT) over transferred data. This allows ping-pong buffers and CRC generation without the CPU.
- **PIO [spec]**: 2 blocks × 4 state machines running deterministic programs clocked from the system clock, with TX
  FIFOs fed by DMA. **PIO is digital-only**: it cannot digitise the AD8232's analog output. In C the ADC still does the
  sampling and PIO only serialises the digital samples.
- **USB CDC on the RP2040**: a full-speed (12 Mbit/s) device with bulk endpoints of up to 64 bytes per packet. The host
  schedules transfers in 1 ms frames (USB 2.0 §5.8). USB bulk has CRC and retry at the link layer, and NAK-based flow
  control. With the Pico SDK + TinyUSB, `tud_task()` has to be called from the main loop. `stdio_usb` writes block up
  to a timeout and then drop.
- **AD8232 [spec]**: a single-supply ECG front end (instrumentation amplifier plus gain/filter stages) whose output is
  referenced to a mid-supply REF. When supplied from 3.3 V its output stays within the ADC's 0–3.3 V range.

## 3. Architecture descriptions (as modelled)

| | A — conventional | B — ADC FIFO + DMA + USB CDC | C — ADC FIFO + DMA + PIO/DMA link |
|---|---|---|---|
| Sample trigger | `time_us_64()` polled in a superloop; `adc_read()` | ADC free-running from CLKDIV | identical to B |
| Sample transfer | CPU | FIFO → DMA ping-pong (2 × 50 samples) | identical to B |
| Packetisation | CPU, after every 50 slots | half-buffer IRQ; packet queued for USB | half-buffer IRQ writes header; DMA chains the block to the PIO TX FIFO |
| Transmission | blocking CDC write (500 ms timeout, then drop) | TinyUSB CDC, serviced by `tud_task()` in the superloop | PIO serial, 2 Mbit/s, 10 bits per byte, unidirectional |
| Flow control | USB NAK → writer blocks → **acquisition stalls** | USB NAK → device buffer fills → **DMA/TX buffer overflow** | none → **receiver buffer overflow** |
| Line errors | corrected by USB retry | corrected by USB retry | detected by CRC → packet discarded |

All architectures send the same packet: magic, sequence number of the first slot (u32), count, flags, 50 × u16
(12-bit code or 0xFFFF for "slot not acquired"), and CRC-16/CCITT-FALSE. That is 110 bytes per 100 ms.

## 4. Mathematical models

**ECG.** Beat *i* has R-time *τᵢ*. RR intervals follow an AR(1) process,
RRᵢ = RR̄(1 + 0.03 zᵢ) with zᵢ = 0.7 zᵢ₋₁ + √0.51 εᵢ.

ECG(t) = Σᵢ Σ_{w∈{P,Q,R,S,T}} a_w exp(−(t − τᵢ − θ_w(RRᵢ))² / 2σ_w²) + W(t)

- θ_P and θ_T scale with √RR.
- W(t) is baseline wander: 0.15 mV at 0.25 Hz plus 0.05 mV at 0.05 Hz.
- Measurement noise is white Gaussian (σ = 0.02 mV) plus 0.01 mV of 50 Hz mains.

The model is continuous in time, so each architecture samples it at its *own* true instants.

**Front end and ADC.**
- AFE output: V = clip(1.65 V + 1 V/mV · ECG, 0, 3.3).
- ADC code: c = clip(⌊V/3.3 · 4096 + n_adc⌋, 0, 4095), with n_adc ~ N(0, 1.5 LSB).
- Reconstruction: V̂ = (c + ½) · 3.3/4096.

**A: CPU timing.**
- Background tasks arrive as a Poisson process at rate λ = ρ / E[D]. Durations D are lognormal and the tasks are
  non-preemptive. Busy intervals come from the vectorised Lindley recursion
  e_k = S_k + max_{j≤k}(a_j − S_{j−1}).
- The loop checks the deadline between tasks. For deadline t_n, the sampling instant is
  s_n = ISR(F(t_n)) + U(0, 0.5 µs) + U(0, 1 µs), where F(t) is the end of the task running at t and ISR(·) adds
  Poisson ISR preemption (2 kHz, 2–10 µs).
- Deadline n is skipped if F(t_n) ≥ t_{n+1}. While the USB write blocks until time u, every deadline whose successor
  also elapsed before u is skipped, and the latest elapsed one is sampled at u.

**B/C: hardware timing.**
- Sampling instant: s_n = nT + j_n, with j_n ~ N(0, 1 ns).
- DMA read instant: r_n = max(r_{n−1}, s_n + 2 µs + δ_n), with δ_n ~ Exp(0.2 µs).
- FIFO overflow at sample n if r_{n−4} > s_n + 2 µs.

**B handoff to USB.** A packet becomes ready at the DMA half-complete time r. It is handed over at
F(r + 5 µs) + 20 µs, i.e. at the next `tud_task()` after the current background task finishes.

**C handoff to PIO.** h = r + 5 µs (IRQ) + 2 µs (DMA chain).

**Transport (per packet, sequential event model).**

*Device buffer* (Q_dev packets):
- On arrival, packets that have already started transmission are released.
- If the buffer is still full, A blocks (with a timeout) and B/C drop the packet.

*USB link:*
- Start time = max(handoff, previous end, host-buffer flow-control release) + U(0, 1 ms) frame alignment.
- The start is pushed to the end of any host polling stall. Stalls arrive as a Poisson process with lognormal
  durations.
- Transmission time = bytes · 8 · 1.25 / 12 Mbit/s.

*PIO link:*
- Start time = max(handoff, previous end).
- Transmission time = bytes · 10 / 2 Mbit/s.
- Arrival at the receiver interface adds 5 µs plus |N(0, 0.5 µs)|.

*Receiver (both links):*
- Arrival = link arrival + host-stack latency (lognormal, median 1 ms), kept in FIFO order.
- The receiver buffer holds Q_rx packets. USB applies back-pressure when it is full; PIO drops.
- Application server: processing time is lognormal (median 0.2 ms). Application stalls arrive as a Poisson process
  with lognormal durations.

*Faults injected identically for all architectures:* packet erasure, corruption (1–3 bit flips), duplication and
pairwise reordering. Physical-line BER applies only to the PIO link (§3).

**Receiver reconstruction.** Each packet's CRC is checked and corrupt packets are discarded. Samples are placed by
sequence number. A slot already seen counts as a duplicate; a slot below the highest sequence number seen counts as
out of order; 0xFFFF counts as missed at acquisition.

**Metrics.**
- Interval error: e_n = (s_n − s_{n'}) − (n − n′)T over consecutive *acquired* samples, using simulator ground-truth
  instants and never receiver timestamps.
- Transport latency: delivery to the receiving application minus packet-ready time.
- Link latency: arrival at the receiver interface minus packet-ready time.
- RMSE vs AFE: received samples compared with the AFE output an ideal sampler would have seen (noise included). This
  isolates quantisation, ADC noise and timing error.
- RMSE filled: gaps are linearly interpolated and compared with the noise-free ECG + wander.
- SNR = 10 log₁₀(var(ref) / MSE).
- R-peaks: an offline Pan-Tompkins-style detector (5–15 Hz band-pass, derivative, squaring, 150 ms integration,
  adaptive threshold, parabolic refinement). Matched to the known R-times within 150 ms (AAMI EC57-style).
- CPU utilisation is **[derived]** from the assumed per-operation costs: A = 100 % (busy-wait), B ≈ 0.025 %,
  C ≈ 0.005 %.

## 5. Experimental methodology

- **Conditions** (overrides in `config.json`):

  | Condition | What changes |
  |---|---|
  | `baseline_no_faults` | no stalls, no BER, no faults; CPU load as in normal |
  | `normal` | ρ = 0.25, tasks with 0.3 ms mean; rare short USB/host stalls; PIO BER 10⁻⁹ |
  | `high_cpu_load` | ρ = 0.75, tasks with 1.5 ms mean, σ = 1.0 |
  | `increased_latency` | host-stack median 15 ms, σ = 0.8 |
  | `reduced_buffer` | Q_dev = 2, Q_rx = 4, plus moderate stalls |
  | `bursty_delays` | USB stalls 0.5/s with 80 ms median; application stalls 0.3/s with 150 ms median |
  | `packet_loss` | 1 % erasure, 0.1 % corruption, 0.1 % duplication, 0.1 % reorder, PIO BER 10⁻⁶ |
  | `receiver_delay` | application processing median 60 ms, stalls 0.1/s with 300 ms median |

- **Trials:** 30 per condition, each a 600 s record at 500 S/s (300 000 samples, 6 000 packets).
- **Seeds:** `SeedSequence([20261016, trial]).spawn(5)` gives the ECG, noise, CPU, hardware and transport streams.
  The *same* trial seeds are used for every architecture and condition (common random numbers). Comparisons are
  therefore **paired**: within a trial, A, B and C see the identical ECG, noise, task workload and host behaviour.
- **Statistics:**
  - Mean ± SD and t-based 95 % CIs per metric.
  - Paired two-sided Wilcoxon signed-rank tests for A–B, B–C and A–C on 14 metrics × 8 conditions (336 tests).
  - Holm–Bonferroni correction across all tests; effect size reported as the rank-biserial correlation.
- **Sensitivity sweeps:** 10 trials × 120 s per point, over CPU load, device buffer size, stall duration and PIO BER.

## 6. Validation

`tests/test_simulation.py` has 16 tests, all passing. They cover:
- expected sample and packet counts
- ADC code limits, the ½-LSB reconstruction bound and AFE clipping
- strictly increasing acquisition instants, contiguous sequence numbers and monotone receiver timestamps
- detection of deliberately dropped packets, 0xFFFF acquisition markers, duplicates and reordering
- CRC detection of all 1–3-bit errors
- device-buffer overflow, receiver overflow without flow control, and back-pressure with USB
- ADC FIFO overflow under forced DMA starvation
- ADC-path RMSE matching the analytic value LSB·√(1/12 + σ²_adc)
- recovery of an injected 5 % loss rate from sequence gaps
- the R-peak detector and scorer on a clean ECG
- bit-exact reproducibility for a fixed seed
- the Holm correction

Checks of consistency with the configuration in the full results:
- **Hardware jitter:** B/C interval-error RMS = 1.413 ns, against a configured 1 ns jitter. The difference of two
  independent N(0, 1 ns) draws should give √2 ns ✓.
- **Packet loss accounting:** in `packet_loss`, 1 998 lost packets were detected for A/B over 30 trials. The
  configured 1 % erasure gives about 1 800, and the 182 CRC-rejected corrupt packets make up the rest ✓.
- **Duplicates and reordering:** 8 450 duplicate *samples* = 169 duplicated packets × 50, matching the configured
  0.1 % ✓.
- **CRC:** every corrupted packet was detected (A/B 182/182, C 312/312).

## 7. Simulation assumptions (and which are specifications)

| Parameter | Value | Status |
|---|---|---|
| ADC resolution, Vref range, conversion time, FIFO depth, CLKDIV formula | 12 bit, 0–3.3 V, 2 µs, 4 | [spec] |
| USB FS bit rate, frame, bulk packet size | 12 Mbit/s, 1 ms, 64 B | [spec] |
| DMA sniffer CRC-16-CCITT | available | [spec] |
| ADC noise | 1.5 LSB RMS Gaussian, no DNL model | [assumption] (meant to stand in for the datasheet's below-12-bit ENOB and erratum E11) |
| ADC sample-clock jitter (B/C) | 1 ns RMS | [assumption] |
| Crystal frequency error | ignored: common to all architectures | [assumption] |
| AFE gain/offset/noise | 1 V/mV, 1.65 V, 20 µV + 10 µV mains | [assumption] |
| CPU background workload | Poisson/lognormal, ρ = 0.25 (0.75 high) | [assumption] |
| Superloop polling granularity, ISRs | 0.5 µs + 1 µs timer resolution; 2 kHz × 2–10 µs | [assumption] |
| Firmware costs (IRQ, tud_task, format) | 5 µs, 20 µs, 30 µs | [assumption] |
| stdio-style write timeout (A) | 500 ms | [assumption]; it resembles the Pico SDK default, but should be checked against the SDK version used |
| USB host-stack latency | lognormal, median 1 ms, σ 0.5 | [assumption] |
| USB/host stall processes | Poisson onset, lognormal duration | [assumption] |
| PIO link | 2 Mbit/s, 10 bits per byte, no flow control, BER 10⁻⁹ | [assumption] (design choice) |
| Receiver path for C | same host-latency/application model as USB; receiver buffer Q_rx = 64 packets | [assumption] |
| Packet size | 50 samples (100 ms) | [assumption] (design choice) |

Simplifications:
- Tasks do not queue behind the sampling code.
- USB link-layer retries add no delay.
- The receiver processes duplicated packets twice.
- Reordering swaps adjacent packets only.
- The ECG is a synthetic morphology. No local ECG dataset was available, so the optional secondary validation on
  recorded ECG was **not performed**.

## 8. Results

All values are mean over 30 trials; [·] is the 95 % CI. Full tables are in `results/summary_tables.md`.

### 8.1 Acquisition timing (figures `01`, `02`, `13`)

| Condition | A interval-error RMS | A max \|error\| | B = C interval-error RMS |
|---|---|---|---|
| baseline / normal | 259.8 µs [259.4, 260.3] | 1 998 µs | 0.0014 µs |
| high_cpu_load | 722.8 µs [722.4, 723.3] | 2 597 µs | 0.0014 µs |

In A, the error is set by how long the background tasks run. It does not change with any transport condition except
when the blocking write stalls the loop: `bursty_delays` and `reduced_buffer` added 0.17 and 0.15 percentage points of
missed samples. B and C were identical, sample for sample, in every condition. The PIO stage has no effect on
acquisition timing, as expected.

### 8.2 Sample loss and buffer overflow (figures `03`–`06`)

| Condition | Acq. loss A | Rx missing A | Rx missing B | Rx missing C | Device overflow B (pkts/10 min) | Receiver overflow C |
|---|---|---|---|---|---|---|
| normal | 0.145 % | 0.145 % | 0 | 0 | 0 | 0 |
| high_cpu_load | **22.26 %** | 22.26 % | 0 | 0 | 0 | 0 |
| reduced_buffer | 0.298 % | 0.303 % | 0.222 % | 0.232 % | 13.3 | 13.9 |
| bursty_delays | 0.318 % | 0.328 % | 0.198 % | **0.007 %** | 11.9 | 0.4 |
| packet_loss | 0.145 % | 1.253 % | **1.110 %** | 1.182 % | 0 | 0 |
| receiver_delay | 0.145 % | 0.145 % | 0 | 0.001 % | 0 | 0.07 |

ADC FIFO overflows were 0 in every trial of B and C. With the assumed DMA latency (0.2 µs mean) and 2 ms between
samples, this is expected rather than evidence.

### 8.3 Transport latency (figures `07`, `08`, `14`, `15`)

| Condition | p99 A | p99 B | p99 C | link-latency jitter SD, B | C |
|---|---|---|---|---|---|
| baseline | 4.08 ms | 4.30 ms | **3.99 ms** | 0.360 ms | 0.0003 ms |
| normal | 4.12 | 4.33 | **4.01** | 0.392 | 0.0003 |
| high_cpu_load | 4.13 | **15.71** | 4.01 | 2.983 | 0.0003 |
| increased_latency | 98.1 | 98.2 | 98.0 (n.s. vs B) | 0.392 | 0.0003 |
| reduced_buffer | 187 | 185 | **157** | 17.8 | 0.0003 |
| bursty_delays | 804 | 803 | **694** | 69.6 | 0.0003 |
| receiver_delay | 1 758 | 1 758 | 1 758 (n.s.) | 0.430 | 0.0003 |

Whenever the host side dominates (`increased_latency`, `receiver_delay`), all three architectures converge.

### 8.4 Waveform fidelity and R-peaks (figures `00`, `09`–`12`)

- **RMSE vs ideally sampled AFE output** isolates acquisition effects:
  - A: 3.01 mV normal, 10.48 mV high load.
  - B = C: 1.230 mV, which equals the analytic quantisation-plus-ADC-noise floor.
  - The extra error in A is the ECG slope × sampling-time error, and it is concentrated on the QRS.
- **RMSE against the noise-free ECG** is dominated by the configured 20 µV × 1000 measurement noise (≈ 20 mV):
  A 21.42, B = C 21.26 mV. It rises only where data is missing, e.g. `packet_loss` 30.1 (B) vs 30.6 (C) mV.
- **R-peak detection:**
  - Sensitivity was ≥ 0.990 and PPV = 1.000 everywhere.
  - Timing-error RMS was 0.772 ms for B/C. This is the detector's own floor at 2 ms sampling, not an acquisition
    effect.
  - For A it was 0.790 ms normally and 1.607 ms under high load.

### 8.5 Sensitivity (figures `16`–`19`)

- **CPU load** (task mean 1 ms, σ = 1):
  - A's acquisition loss grows from 1.9 % (ρ = 0.1) to 16.2 % (ρ = 0.85), and its RMS interval error from 319 to
    710 µs.
  - B's p99 latency grows from 5.6 to 12.0 ms.
  - C does not change (3.98 ms).
- **Device buffer** under bursty stalls:
  - B's missing rate falls from 3.2 % (1 packet) to 0.02 % (16) and 0 (32).
  - C stays at 0.21 %, set by its *receiver* buffer (8 packets in this sweep), not the device buffer.
  - **B beats C once Q_dev ≥ 16**, because USB flow control turns stalls into delay rather than loss.
- **Stall duration** (USB and application stalls at 0.5/s):
  - The A/B missing rate reaches 6.3 % at a 0.4 s median.
  - C loses nothing, because its 64-packet receiver buffer (6.4 s) absorbs the stalls.
- **PIO BER:**
  - C's losses are negligible at ≤ 10⁻⁷ but reach 0.1 % at 10⁻⁶, 0.8 % at 10⁻⁵ and 8.6 % at 10⁻⁴.
  - A and B are unaffected, because of USB retry.

## 9. Statistical analysis

All B-vs-C acquisition metrics (interval error, acquisition loss, RMSE vs AFE) have identical paired values,
so p = 1. **Significant B-vs-C differences** (Holm-adjusted p < 0.05; full list in `paired_tests.csv` /
`summary_tables.md`):

- **Favouring C:**
  - Mean, p99 and jitter of latency, and link-latency jitter, in baseline, normal, high_cpu, reduced_buffer,
    bursty and packet_loss (p_Holm ≈ 6·10⁻⁷). The size of the improvement depends on the condition:
    - normal: 0.32 ms at p99, i.e. 7 %
    - high_cpu: 11.7 ms
    - bursty: 110 ms
  - Missing samples and device overflow under bursty stalls (p_Holm = 2·10⁻³).
  - RMSE-filled under bursty stalls; R-peak sensitivity under bursty stalls (p_Holm = 0.0496, borderline).
- **Favouring B:**
  - Missing samples and RMSE-filled under `packet_loss` (BER 10⁻⁶), p_Holm = 3·10⁻⁴.
  - Receiver overflow under reduced_buffer.
- **No significant difference:**
  - total missing samples under reduced_buffer and receiver_delay
  - p99 latency under increased_latency and receiver_delay
  - R-peak timing error everywhere

A vs B: B is significantly better (p_Holm ≤ 3·10⁻⁴) on interval error, acquisition loss and RMSE vs AFE in every
condition. B's p99 latency is *worse* than A's in normal (+0.2 ms) and high_cpu (+11.6 ms), because in this model B
hands packets to `tud_task()` in the superloop while A writes straight after sampling. Remember that A's own latency
is measured from a sampling instant that may already be late.

## 10. Discussion — best architecture per criterion

| Criterion | Winner (this model) | Basis |
|---|---|---|
| Analog ECG acquisition timing | **B = C (tie)** | ns vs hundreds of µs; B and C use the identical ADC path |
| End-to-end data integrity | **depends on the fault model; no single winner** | C under host bursts with a large receiver buffer; B under line bit errors or with a large device buffer; equal with small buffers |
| Low-latency streaming | **C** | lowest p50/p99 and ~1000× lower link jitter in all non-host-dominated conditions; differences ≤ 0.3 ms in normal conditions |
| Robustness under load | **C** (B close behind) | CPU load does not affect C at all; B keeps acquisition intact but its latency grows; A loses up to 22 % of samples |
| Overall suitability for ECG on RP2040 | **B or C over A; B is the pragmatic default** | see below |

The decisive result is that **hardware-paced acquisition (B or C) is essential and A is not suitable** under any
non-trivial CPU load. Between B and C:

- For a 500 S/s ECG, the latency gains from PIO (≈ 0.3 ms at p99 normally) are clinically irrelevant for display or
  monitoring.
- C's real advantage is *decoupling from the CPU and USB stack*. That matters if the firmware also runs heavy
  processing, or if hard real-time delivery to dedicated hardware is required.
- C needs a second receiver device, has no flow control or retransmission (so link noise turns into data loss), and
  in this model its integrity advantage under bursts depends on giving the receiver a larger buffer than B's device
  buffer.
- B can match C's integrity simply by increasing its device buffer to ≥ 16 packets (about 1.8 kB of RAM).
- The B-vs-C integrity comparison is therefore largely a comparison of *buffer placement and flow-control policy*,
  not of PIO itself.

## 11. Limitations and what hardware must establish

1. **Assumed distributions decide the transport results.** These include USB host latency and stalls, firmware
   costs, CPU workload and BER. The rankings under bursty, reduced_buffer and packet_loss reverse when the
   buffer/BER assumptions change (§8.5). A defensible winner for *integrity* needs:
   (a) measured host-side CDC delivery-latency and stall distributions on the target PC/OS;
   (b) the measured BER of the intended PIO cable or link;
   (c) the real application CPU workload profile.
2. The ADC is modelled ideally plus Gaussian noise. Real ENOB, DNL (E11) and AD8232 noise set the waveform-error
   floor for all architectures equally. They must be measured but are unlikely to change the ranking.
3. B/C sampling jitter (1 ns) is an assumption. It would need measuring with a logic analyser on an ADC-synchronous
   GPIO, but it is orders of magnitude below what matters for ECG.
4. USB retries, hub contention, CDC driver buffering and TinyUSB internals are abstracted into a single
   stall + latency model.
5. The ECG is synthetic: no arrhythmias, motion artefacts or electrode-off events, and no real-record validation.
6. CPU utilisation figures are derived from assumed costs, not measured.

`HARDWARE_ROADMAP.md` describes the measurements, including firmware based on `pico-examples/adc/dma_capture`.

## 12. Conclusion

Under identical ECG input, sampling rate, quantisation and packet format, this simulation finds:

1. **Hardware-paced ADC + DMA acquisition makes a decisive difference.** RMS interval error falls by about five
   orders of magnitude, from 260–723 µs to about 1 ns, and acquisition loss falls from up to 22 % to 0.
2. **PIO + DMA transmission (C vs B) gives a statistically measurable improvement in transmission-timing
   determinism.** Link-latency jitter drops from 0.36–70 ms to 0.3 µs, and latency becomes immune to CPU load.
3. **It does not give a general improvement in data integrity.** C is better under host stall bursts, worse under
   link bit errors, and equal with small buffers. All of these depend on buffer and flow-control choices that B can
   also make.
4. **Acquisition-level ECG quality is identical for B and C.**

**Final answer:** Yes for timing and latency of the transmission stage. No demonstrable, assumption-independent
improvement for ECG data integrity or acquisition timing. Confirming either result requires the hardware
measurements in §11.

## References

1. Raspberry Pi Ltd, *RP2040 Datasheet* (ADC §4.9, DMA §2.5, PIO §3, USB §4.1, errata incl. RP2040-E11). https://datasheets.raspberrypi.com/rp2040/rp2040-datasheet.pdf
2. Raspberry Pi Ltd, *Raspberry Pi Pico C/C++ SDK* (hardware_adc, hardware_dma, hardware_pio, pico_stdio_usb). https://datasheets.raspberrypi.com/pico/raspberry-pi-pico-c-sdk.pdf
3. Raspberry Pi Ltd, pico-examples `adc/dma_capture/dma_capture.c`. https://github.com/raspberrypi/pico-examples/blob/master/adc/dma_capture/dma_capture.c
4. Analog Devices, *AD8232 Single-Lead, Heart Rate Monitor Analog Front End*, datasheet. https://www.analog.com/en/products/ad8232.html
5. USB Implementers Forum, *Universal Serial Bus Specification Revision 2.0* (2000), §5.8 bulk transfers, §8 protocol layer.
6. USB-IF, *Class Definitions for Communications Devices (CDC) 1.2* / PSTN subclass.
7. TinyUSB project documentation. https://docs.tinyusb.org
8. J. Pan, W. J. Tompkins, "A Real-Time QRS Detection Algorithm," *IEEE Trans. Biomed. Eng.*, BME-32(3), 230–236, 1985.
9. P. E. McSharry, G. D. Clifford, L. Tarassenko, L. A. Smith, "A dynamical model for generating synthetic electrocardiogram signals," *IEEE Trans. Biomed. Eng.*, 50(3), 289–294, 2003.
10. ANSI/AAMI EC57, *Testing and reporting performance results of cardiac rhythm and ST segment measurement algorithms*.
11. S. Holm, "A simple sequentially rejective multiple test procedure," *Scand. J. Statist.*, 6(2), 65–70, 1979.
12. F. Wilcoxon, "Individual comparisons by ranking methods," *Biometrics Bulletin*, 1(6), 80–83, 1945.
