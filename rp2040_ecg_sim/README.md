# RP2040 + AD8232 ECG acquisition — simulation study

Monte Carlo, packet-level discrete-event simulation comparing three acquisition/transmission
architectures for an AD8232 ECG front end on a Raspberry Pi Pico (RP2040):

| | Acquisition | Transmission |
|---|---|---|
| **A** | CPU superloop polls a timer and calls `adc_read()` | blocking USB CDC write (stdio-style, 500 ms timeout) |
| **B** | ADC free-running (CLKDIV) → 4-deep FIFO → DMA ping-pong | USB CDC (TinyUSB task in main loop) |
| **C** | identical to B | DMA → PIO state machine → custom unidirectional serial link |

> **All numbers produced by this project are simulation outputs.** No hardware was measured.
> See `REPORT.md` §7 for which parameters are documented hardware specifications and which are assumptions.

## Setup

Python ≥ 3.10.

```bash
pip install numpy scipy pandas matplotlib pytest
```

Versions used for the committed results are in `results/manifest.json`.

## Run

```bash
python -m pytest -q tests          # 16 validation tests (~3 s)
python run_simulation.py --quick   # smoke run: 3 trials x 60 s -> results_quick/
python run_simulation.py           # full study: 8 conditions x 30 trials x 10 min + sensitivity sweeps -> results/
```

The full run takes a few minutes on 4 cores (`n_workers` in `config.json`).

## Layout

```
config.json              all parameters, conditions and sensitivity sweeps
run_simulation.py        main script (experiment, statistics, plots, manifest)
src/ecgsim/
  ecg.py                 synthetic ECG (P,QRS,T, HRV, wander, noise), AD8232 + 12-bit ADC model
  timing.py              Poisson busy/stall interval primitives (vectorised Lindley recursion)
  architectures.py       A/B/C acquisition, DMA/FIFO, packetisation, fault injection
  pipeline.py            device TX buffer -> USB CDC or PIO link -> receiver buffer -> receiver app
  packets.py             common binary packet format + CRC-16/CCITT
  receiver.py            sequence-number reconstruction, gap/dup/out-of-order/CRC detection
  rpeak.py               Pan-Tompkins-style detector + EC57-style scoring
  analysis.py            aggregation, 95% CIs, paired Wilcoxon tests + Holm correction
  plots.py               figures
tests/test_simulation.py validation tests
results/csv/             trial_metrics.csv, aggregate_metrics.csv, paired_tests.csv, sensitivity_*.csv
results/plots/           PNG figures (200 dpi)
results/manifest.json    reproducibility manifest (versions, seeds, date)
REPORT.md                research report
HARDWARE_ROADMAP.md      plan for validating on a real Pico + AD8232
```

## Interpreting the results

* `interval_err_*` — error of the **true acquisition instants** (ground truth inside the simulator),
  never the host arrival time. Host delivery variability shows up only in `lat_*`.
* `acq_sample_loss_rate` — samples never acquired (A: missed deadlines; B/C: ADC FIFO overflow).
  `rx_missing_rate` — samples absent at the receiver for any reason (acquisition, buffer overflow, link loss, CRC reject).
* `lat_*` — packet-ready → delivered to receiver application. `link_lat_*` — packet-ready → arrival at the receiver
  interface (isolates the device + link stage, which is where PIO can make a difference).
* `rmse_vs_afe_mv` — received samples vs what an ideal sampler would have read from the AFE (isolates quantisation,
  ADC noise and sampling-time error). `rmse_filled_mv` — gaps linearly interpolated, vs the noise-free ECG.
* `cpu_util_est` — derived from assumed per-operation costs; not a measurement.
* `paired_tests.csv` — paired by trial (same ECG and acquisition random draws for all architectures). A difference is
  called only when Holm-adjusted p < 0.05.

Changing an assumption: edit `config.json` (or add a condition) and rerun. Every condition is a set of overrides
on the base parameters, e.g. `"high_cpu_load": {"cpu.load": 0.75, ...}`.
