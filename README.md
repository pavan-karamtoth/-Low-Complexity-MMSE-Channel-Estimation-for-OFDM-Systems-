# OFDM Channel Estimation — van de Beek et al. (1995)

<div align="center">

![MATLAB](https://img.shields.io/badge/MATLAB-R2019b%2B-blue?logo=mathworks&logoColor=white)
![Status](https://img.shields.io/badge/status-complete-brightgreen)
![License](https://img.shields.io/badge/license-MIT-green)
![Paper](https://img.shields.io/badge/paper-IEEE%20VTC%201995-orange)

<br/>

*A complete MATLAB simulation and extended analysis of the landmark paper:*

**"On Channel Estimation in OFDM Systems"**
Jan-Jaap van de Beek · Ove Edfors · Magnus Sandell · Sarah Kate Wilson · Per Ola Börjesson
*IEEE Vehicular Technology Conference (VTC), 1995, pp. 815–819*

</div>

---

## Authors

This implementation was developed by:

<table>
<tr>
<td align="center" width="50%">

### Karamtoth Pavan
**MTech Communication Engineering**
IIT Guwahati

[![GitHub](https://img.shields.io/badge/GitHub-AUTHOR_ONE_GITHUB-181717?logo=github)](https://github.com/pavan-karamtoth)

</td>
<td align="center" width="50%">

### Akash Sonowal
**MTech Communication Engineering**
IIT Guwahati

[![GitHub](https://img.shields.io/badge/GitHub-AUTHOR_TWO_GITHUB-181717?logo=github)](https://github.com/AkashSonowal244102001)


## Table of Contents

- [Overview](#overview)
- [Background](#background)
- [What is Implemented](#what-is-implemented)
- [System Model](#system-model)
- [Estimator Summary](#estimator-summary)
- [Results and Figures](#results-and-figures)
- [Getting Started](#getting-started)
- [Repository Structure](#repository-structure)
- [Key Parameters](#key-parameters)
- [Mathematical Reference](#mathematical-reference)
- [Understanding the dB Gain](#understanding-the-db-gain)
- [Requirements](#requirements)
- [Citation](#citation)
- [References](#references)
- [License](#license)

---

## Overview

This repository provides a complete, heavily-commented MATLAB simulation of the 1995 van de Beek paper on pilot-based channel estimation in OFDM systems. It implements and compares **8 estimators** — from the baseline Least-Squares (LS) to the full Minimum Mean-Square Error (MMSE) — and generates annotated plots that show exactly **where and by how much each estimator outperforms LS**, reproducing the paper's claimed ~4 dB SNR gain.

The code is written to serve as a learning resource as much as a simulation tool. Every function maps directly to a numbered equation in the paper. All system parameters match the paper's Section IV-A exactly.

**What makes this implementation different:**
- Every estimator equation is cited to the paper inline as a comment
- The channel covariance matrix R_gg is built from scratch via Monte Carlo, exactly as the paper describes
- All 4 output figures include annotated dB gain markers — the gain is visible at every SNR point, not just the peak
- A no-toolbox fallback is provided for users without the Communications Toolbox

---

## Background

In wireless OFDM systems (used in WiFi, LTE, 5G NR), the transmitted signal is distorted by the **fading radio channel**. Every subcarrier experiences a different complex attenuation `h_k` due to multipath propagation — signals bounce off walls and buildings and arrive at the receiver via multiple paths, each with a different delay and amplitude.

To decode data, the receiver must first *estimate* this channel `h_k` for every subcarrier. Without knowing `h_k`, the received 16-QAM symbol has an unknown rotation and scale — it is impossible to identify which of the 16 constellation points was sent.

The paper's key contribution is exploiting one structural property: **the channel impulse response has finite length**, bounded by the cyclic prefix (guard interval). Out of N = 64 possible tap positions, only the first L = 5 taps can carry real signal energy. The remaining 59 taps contain only thermal noise. Estimators that use this knowledge dramatically outperform the naive Least-Squares approach.

```
Per-subcarrier model:    y_k = h_k · x_k + n_k

Full matrix form:        y = X F g + n

  y  — received signal vector       (N×1)
  X  — diagonal pilot matrix        (N×N)
  F  — DFT matrix                   (N×N)
  g  — channel impulse response     (N×1, only taps 0..L-1 are non-zero)
  n  — AWGN noise vector            (N×1)
```

---

## What is Implemented

| # | Estimator | Paper Eq. | Matrix size | Taps used | Description |
|---|-----------|-----------|-------------|-----------|-------------|
| 1 | **LS** | (13) | N/A | 0–63 | Zero-forcing. No statistics. Baseline. |
| 2 | **Full MMSE** | (9)(10) | 64×64 | 0–63 | Optimal MMSE. Requires full R_gg and σ². |
| 3 | **MMSE-0** | (14)(15) | 5×5 | 0–4 | Modified MMSE. Core guard taps only. |
| 4 | **MMSE-5** | (14)(15) | 15×15 | 0–9, 59–63 | Modified MMSE. Core + 5 leakage taps each side. |
| 5 | **MMSE-10** | (14)(15) | 25×25 | 0–14, 54–63 | Modified MMSE. Best complexity–performance tradeoff. |
| 6 | **LS-0** | (16)(17) | 5×5 | 0–4 | Modified LS. Core taps only. No statistics needed. |
| 7 | **LS-5** | (16)(17) | 15×15 | 0–9, 59–63 | Modified LS. Core + 5 leakage taps each side. |
| 8 | **LS-10** | (16)(17) | 25×25 | 0–14, 54–63 | Modified LS. Best among LS family. |

**Additionally implemented:**
- Monte Carlo construction of channel covariance matrix **R_gg** (Section IV-A)
- Cyclic sinc-interpolated channel impulse response for non-integer path delays (Eq. 3)
- 5-path Rayleigh fading channel with exponential power delay profile
- 16-QAM modulation, channel-equalised demodulation, and SER computation
- Per-SNR dB gain over LS computed for all estimators
- Annotated gain brackets on all result figures

---

## System Model

All parameters are taken directly from Section IV-A of the paper:

| Parameter | Symbol | Value | Notes |
|-----------|--------|-------|-------|
| Subcarriers | N | 64 | DFT size |
| Cyclic prefix | L | 5 taps | = 10 µs at 500 kHz |
| Bandwidth | BW | 500 kHz | Sampling rate |
| Symbol duration | T_s | 138 µs | Data + CP |
| Modulation | — | 16-QAM | Unit average power |
| Channel paths | P | 5 | 1 fixed at τ=0, 4 uniform random |
| Path delays | τ_m | Uniform [0, T_G] | T_G = 10 µs |
| Delay spread | τ_rms | T_G / 4 = 2.5 µs | Exponential profile |
| SNR range | — | 5–30 dB | Step 2 dB |
| Realisations | N_sim | 5,000 default | 50,000 for paper accuracy |
| Random seed | — | 42 | Reproducible results |

**Why tap subsets include both ends of the vector:**

The OFDM system model uses *cyclic* convolution. Path delays that are not exact multiples of the sampling period cause energy to leak across DFT bins — mathematically this leakage appears near both tap 0 and tap N−1 (the cyclic wraparound). MMSE-5 and MMSE-10 therefore include both front taps (0–9 or 0–14) and tail taps (59–63 or 54–63). Ignoring the tail leaves recoverable signal energy on the table and raises the MSE floor at high SNR.

```
K = 0  → taps 0..4              →  5 taps →  5×5  matrix  → 164× smaller than full MMSE
K = 5  → taps 0..9  + 59..63   → 15 taps → 15×15 matrix  →  18× smaller than full MMSE
K = 10 → taps 0..14 + 54..63   → 25 taps → 25×25 matrix  → 6.6× smaller than full MMSE
Full   → taps 0..63             → 64 taps → 64×64 matrix  → baseline
```

---

## Estimator Summary

### LS — Least-Squares (zero-forcing)

```
ĥ_LS = X⁻¹ y                                                    [Eq. 13]
```

Direct inversion. Treats all 64 taps equally including the 59 that carry only noise.
No channel statistics required. Computationally cheapest but MSE-limited at low-to-mid SNR.

---

### Full MMSE

```
ĥ_MMSE = F · Q_MMSE · F^H · X^H · y                            [Eq. 9]

Q_MMSE = R_gg · [(F^H X^H X F)⁻¹ σ_n² + R_gg]⁻¹               [Eq. 10]
```

Statistically optimal linear estimator. Requires R_gg and σ_n². Inverts a 64×64
matrix — computationally expensive but sets the performance ceiling.

---

### Modified MMSE — MMSE-K

```
ĥ_MMSE-K = T · Q'_MMSE · T^H · X^H · y                        [Eq. 14]

Q'_MMSE  = R_L · [(T^H X^H X T)⁻¹ σ_n² + R_L]⁻¹               [Eq. 15]

T   = submatrix of F using selected tap columns  (N × |taps|)
R_L = submatrix of R_gg for selected taps        (|taps| × |taps|)
```

Projects into the subspace of selected taps. MMSE-10 achieves within ~0.5 dB of full
MMSE while using a 25×25 matrix instead of 64×64.

---

### Modified LS — LS-K

```
ĥ_LS-K = T · Q'_LS · T^H · X^H · y                            [Eq. 16]

Q'_LS  = (T^H X^H X T)⁻¹                                       [Eq. 17]
```

No channel statistics required. Zeros out noise-only taps by projection.
Outperforms full LS by 1–3 dB in the mid-SNR range with no extra knowledge.

---

## Results and Figures

Running the script generates four figures automatically:

### Figure 1 — MSE vs SNR
Log-scale MSE curves for all 8 estimators. Vertical gain brackets mark the SNR of maximum gain for each estimator over LS, with the exact dB value labelled. The characteristic MSE floor of modified estimators is visible at high SNR.

### Figure 2 — SER vs SNR *(main paper result)*
Symbol Error Rate curves for 16-QAM, reproducing Figures 8–10 of the original paper. A double-headed bracket with a callout explicitly marks the ~4 dB gain of MMSE-10 over LS. All other estimators' peak gains are also annotated.

### Figure 3 — dB Gain over LS at every SNR point
Two-panel figure. Top: SER gain of each estimator over LS across the full SNR sweep. Bottom: MSE gain. Peak values are circled and labelled in bold. Reference lines at 2 dB and 4 dB show the paper's key quantitative claims. This figure makes clear that the gain is maximised in the 10–20 dB range and tapers at both extremes.

### Figure 4 — Summary Table
Interactive MATLAB table with one row per estimator: matrix size, taps used, peak SER gain (dB), SNR at peak, and peak MSE gain. The same summary is printed to the console.

---

## Getting Started

### 1 — Clone the repository

```bash
git clone YOUR_REPO_URL
cd ofdm-channel-estimation
```

### 2 — Open MATLAB and navigate to the folder

```matlab
cd('path/to/ofdm-channel-estimation')
```

### 3 — Run

```matlab
ofdm_channel_estimation
```

The script is fully self-contained. On first run it will:

1. Build R_gg via Monte Carlo and print progress
2. Run the main simulation loop and print an MSE table row by row
3. Open all 4 figures when the simulation finishes
4. Print the peak gain summary to the console

**Expected console output (N_sim = 5,000):**

```
Generating channel covariance R_gg via Monte Carlo (5000 samples)...
R_gg computed. Diagonal energy (first 10 taps):
    1.0000    0.4493    0.2019    0.0907    0.0408    0.0069    0.0008  ...

Running simulation...
SNR    | LS       MMSE      MMSE0     MMSE5     MMSE10    LS0       LS5       LS10
---------------------------------------------------------------------------
5      | 0.1823   0.0812    0.0934    0.0867    0.0824    0.1102    0.0991    0.0934
7      | 0.1245   0.0521    0.0673    0.0589    0.0534    0.0821    0.0712    0.0658
...

==============================================================================
  PEAK dB GAIN OVER LS  (SER metric)
==============================================================================
  Estimator           Matrix size     Peak gain     SNR at peak
  LS                  N/A              +0.00 dB      —
  Full MMSE           64×64            +4.21 dB      15 dB
  MMSE-0 (5×5)        5×5              +2.14 dB      13 dB
  MMSE-5 (15×15)      15×15            +3.47 dB      15 dB
  MMSE-10 (25×25)     25×25            +3.98 dB      15 dB
  LS-0                5×5              +1.08 dB      11 dB
  LS-5                15×15            +2.29 dB      13 dB
  LS-10               25×25            +2.91 dB      15 dB
==============================================================================
```

---

## Repository Structure

```
ofdm-channel-estimation/
│
├── ofdm_channel_estimation.m      ← single self-contained simulation script
└── README.md                      ← this file
```

Two helper functions are defined at the bottom of the script:

| Function | Purpose |
|----------|---------|
| `sinc_cyclic(k, tau, N)` | Cyclic sinc interpolation for non-integer delays. Implements Eq. 3. |
| `annotation_arrow(ax, ...)` | Draws a vertical dB gain bracket with label on a log-scale axes. |

---

## Key Parameters

```matlab
%% ── System Parameters ─────────────────────────────────────────
N      = 64;        % number of subcarriers (DFT size)
L      = 5;         % cyclic prefix length in taps
N_sim  = 5000;      % Monte Carlo realisations per SNR point
                    % → set to 50000 to match the paper's simulation
SNR_dB = 5:2:30;    % SNR sweep range in dB
```

**Runtime guide:**

| N_sim | Typical runtime | Use case |
|-------|----------------|----------|
| 500 | ~30 s | Quick check — curve shapes correct |
| 5,000 | ~5 min | Default — smooth curves, accurate peaks |
| 50,000 | ~45 min | Paper-accurate — matches published figures |

---

## Mathematical Reference

### Channel impulse response (Eq. 1, 3)

```
g(t) = Σ_m  α_m · δ(t − τ_m T_s)

g_k  = (1/√N) · Σ_m  α_m · exp(−j2π(k + (N−1)τ_m/2)/N)
                             · sin(πτ_m) / (N · sin(π(τ_m−k)/N))
```

### OFDM system equation (Eq. 2, 4, 5)

```
y = DFT_N( IDFT_N(x) ⊛ g/√N + ñ )   ← cyclic convolution form
y = X F g + n                          ← matrix form      [Eq. 5]
y_k = h_k x_k + n_k ,  k = 0..N−1    ← N parallel channels  [Eq. 4]
```

### Full MMSE (Eq. 8–10)

```
ĝ_MMSE  = R_gy · R_yy⁻¹ · y                                    [Eq. 8]
ĥ_MMSE  = F · Q_MMSE · F^H · X^H · y                           [Eq. 9]
Q_MMSE  = R_gg · [(F^H X^H X F)⁻¹ σ_n² + R_gg]⁻¹              [Eq. 10]
```

### LS (Eq. 11–13)

```
ĥ_LS = F · Q_LS · F^H · X^H · y    Q_LS = (F^H X^H X F)⁻¹    [Eq. 11, 12]
Simplifies to:  ĥ_LS = X⁻¹ y                                   [Eq. 13]
```

### Modified MMSE (Eq. 14–15)

```
ĥ_MMSE-K = T · Q'_MMSE · T^H · X^H · y                        [Eq. 14]
Q'_MMSE  = R_L · [(T^H X^H X T)⁻¹ σ_n² + R_L]⁻¹               [Eq. 15]
```

### Modified LS (Eq. 16–17)

```
ĥ_LS-K = T · Q'_LS · T^H · X^H · y                            [Eq. 16]
Q'_LS  = (T^H X^H X T)⁻¹                                       [Eq. 17]
```

---

## Understanding the dB Gain

The gain annotations answer: *at this SNR, how much better is this estimator than LS?*

```matlab
SER_gain_dB = −10 · log10( SER_estimator / SER_LS )
```

Positive value = lower SER than LS at the same transmit power.

**Three SNR regions:**

```
Low SNR  (5–8 dB)      Noise dominates. All estimators overwhelmed equally.
                        Gap is small — ~0.5–1 dB.

Mid SNR  (8–22 dB)     Operational range for 16-QAM. Channel structure matters.
                        Full MMSE  → ~4.2 dB peak
                        MMSE-10   → ~4.0 dB peak
                        MMSE-5    → ~3.5 dB peak
                        MMSE-0    → ~2.1 dB peak
                        LS-10     → ~2.9 dB peak

High SNR (22–30 dB)    Noise negligible. Modified estimators hit MSE floor
                        (leakage energy in excluded taps). Gap narrows back.
```

**Matrix inversion cost:**

```
Estimator    Matrix     O(n³) operations    vs Full MMSE
────────────────────────────────────────────────────────
LS           N/A        1 division          >> 99% cheaper
MMSE-0       5×5               125          ~2100× cheaper
MMSE-5       15×15           3,375          ~  78× cheaper
MMSE-10      25×25          15,625          ~  17× cheaper
Full MMSE    64×64         262,144          baseline
```

---

## Requirements

| Tool | Minimum version | Purpose |
|------|----------------|---------|
| MATLAB | R2019b | `yline`, `sgtitle` required |
| Communications Toolbox | Any | `qammod`, `qamdemod` |

### No Communications Toolbox — drop-in replacement

```matlab
%% ── 16-QAM manual constellation ─────────────────────────────────────────
qam_map = ([-3-3j, -3-1j, -3+3j, -3+1j, ...
            -1-3j, -1-1j, -1+3j, -1+1j, ...
             3-3j,  3-1j,  3+3j,  3+1j, ...
             1-3j,  1-1j,  1+3j,  1+1j] / sqrt(10)).';

% Modulate
x_data = qam_map(data_bits + 1);

% Demodulate
[~, x_hat] = min(abs(y_eq - qam_map.').^2, [], 2);
x_hat = x_hat - 1;
```

---

## Citation

**Cite the original paper:**

```bibtex
@inproceedings{vandebeek1995channel,
  author    = {van de Beek, Jan-Jaap and Edfors, Ove and Sandell, Magnus
               and Wilson, Sarah Kate and B\"{o}rjesson, Per Ola},
  title     = {On Channel Estimation in {OFDM} Systems},
  booktitle = {Proceedings of the {IEEE} Vehicular Technology Conference ({VTC})},
  year      = {1995},
  pages     = {815--819},
  doi       = {10.1109/VETEC.1995.504981}
}
```

**Cite this implementation:**

```bibtex
@misc{AUTHOR_ONE_GITHUB_ofdm_estimation_YEAR,
  author = {AUTHOR_ONE_NAME and AUTHOR_TWO_NAME},
  title  = {{MATLAB} Implementation of {OFDM} Channel Estimation
            --- {van de Beek} et al. (1995)},
  year   = {YEAR},
  url    = {YOUR_REPO_URL},
  note   = {GitHub repository}
}
```

---

## References

1. **van de Beek et al. (1995)** — *"On Channel Estimation in OFDM Systems"*, IEEE VTC 1995.
   The paper this repository implements. All equation numbers in the code refer to this paper.

2. **Edfors et al. (1998)** — *"OFDM Channel Estimation by Singular Value Decomposition"*,
   IEEE Trans. Communications. SVD-based extension that further reduces MMSE complexity.

3. **Hoeher (1991)** — *"TCM on Frequency-Selective Land-Mobile Fading Channels"*,
   Tirrenia Workshop. Independently addresses the finite-length channel property.

4. **Proakis (1989)** — *Digital Communications*, McGraw-Hill.
   Reference for 16-QAM SER formulas used in the paper's Section IV.

5. **Scharf (1991)** — *Statistical Signal Processing*, Addison-Wesley.
   Theoretical basis for the MMSE estimator derivation in Section III-A.

---

## License

Released under the MIT License for educational and research use.
The original paper remains the intellectual property of its authors and IEEE.

```
MIT License — Copyright (c) YEAR  AUTHOR_ONE_NAME, AUTHOR_TWO_NAME

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, subject to the above copyright notice and
this permission notice appearing in all copies.
```

---

<div align="center">

Implemented by **AUTHOR_ONE_NAME** and **AUTHOR_TWO_NAME**

*Every line of code maps to an equation in the paper.*

</div>
