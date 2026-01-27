# 2×2 MIMO-OFDM Simulation (ZF vs ZF-SIC + MRC)

## Overview
This project simulates a **2×2 MIMO-OFDM** wireless link over a **frequency-selective Rayleigh fading** channel with **AWGN**. It compares two receiver strategies:

1. **ZF-only** (per-subcarrier zero-forcing detection of both spatial streams)
2. **ZF-SIC + MRC** (detect stream-2 using ZF, cancel its interference, then detect stream-1 using MRC)

The script produces **BER vs SNR** curves for **BPSK** or **QPSK**.

---

## Features
- True OFDM processing:
  - Subcarrier mapping
  - **IFFT/FFT**
  - **Cyclic Prefix (CP) insertion/removal**
- Frequency-selective **L-tap Rayleigh** 2×2 MIMO channel
- Perfect timing and perfect CSI (uses true channel frequency response)
- Modulation options:
  - BPSK (1 bit/symbol)
  - QPSK (2 bits/symbol, Gray mapping)
- Output: BER performance comparison for ZF-only vs ZF-SIC+MRC

---

## Parameters (default)
- `Nfft = 64`  : number of subcarriers (FFT size)
- `Ncp  = 16`  : cyclic prefix length (samples)
- `Nsym = 1500`: number of OFDM symbols
- `Lch  = 6`   : channel taps (multipath length)
- `SNR_dB = 0:2:24`

**CP condition (for ISI-free OFDM):**
`Ncp >= Lch - 1`  → with `Ncp=16` and `Lch=6`, the condition is satisfied.

---

## How it Works (High-Level)
### Transmitter
1. Generate random bits for each transmit antenna
2. Map bits to BPSK/QPSK symbols per subcarrier
3. Perform IFFT (OFDM modulation)
4. Add cyclic prefix
5. Serialize per-transmit-antenna waveform

### Channel
- Apply a frequency-selective 2×2 Rayleigh channel (time-domain convolution)
- Add complex AWGN at desired SNR

### Receiver
1. Reshape waveform into OFDM symbols
2. Remove CP and apply FFT
3. Per-subcarrier detection:
   - **ZF-only:** detect both streams directly
   - **ZF-SIC+MRC:** detect stream-2 → cancel → MRC detect stream-1
4. Demap symbols to bits and compute BER

---

## Running the Simulation
1. Open MATLAB and run the main script.
2. Choose modulation when prompted:
   - `1` → BPSK
   - `2` → QPSK
3. The script prints BER results and plots BER vs SNR.

---

## Output
A BER vs SNR plot comparing:
- **ZF-only**
- **ZF-SIC + MRC**

---

## Notes / Assumptions
- Perfect CSI and perfect timing are assumed.
- Results improve (smoother curves) by increasing `Nsym`.
- For more realism, add pilots + channel estimation (LS/MMSE) and timing/CFO impairments.

---
