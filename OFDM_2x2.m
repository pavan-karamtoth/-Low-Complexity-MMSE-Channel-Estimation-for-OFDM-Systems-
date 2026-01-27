%% 2x2 MIMO-OFDM (1 user, spatial multiplexing) with ZF-only vs ZF-SIC(+MRC)
% - True OFDM: subcarriers + IFFT/FFT + cyclic prefix
% - Frequency-selective Rayleigh channel (L-tap), perfect CSI at RX
% - Modulation: BPSK or QPSK (Gray)
%
% Receiver 1: per-subcarrier ZF detect both streams (ZF-only)
% Receiver 2: per-subcarrier ZF -> detect stream-2 -> SIC cancel -> MRC detect stream-1 (ZF-SIC+MRC)
%
% Output: BER vs SNR (for both BPSK and QPSK, always BER over bits)

clc; clear; close all;

%% -------------------- User Inputs --------------------
mod = input("Enter modulation: (1-->BPSK) (2-->QPSK): ");
if mod ~= 1 && mod ~= 2
    error('Choose mod = 1 (BPSK) or 2 (QPSK).');
end

noTx = 2; noRx = 2;

Nfft = 64;        % number of subcarriers
Ncp  = 16;        % cyclic prefix length
Nsym = 1500;      % OFDM symbols per frame (increase for smoother curves)
Lch  = 6;         % channel taps

SNR_dB = 0:2:24;  % SNR per received complex sample (approx)
SNRlin = 10.^(SNR_dB/10);

% Power delay profile (normalized)
pdp = exp(-(0:Lch-1));
pdp = pdp / sum(pdp);

% bits per symbol
if mod == 1
    bps = 1;
else
    bps = 2;
end

%% -------------------- Generate Bits & Map to OFDM Grid --------------------
% Bits per Tx stream
Nb = Nfft*Nsym*bps;

txBits = randi([0 1], Nb, noTx);  % [Nb x 2]

% Map bits -> symbols (unit average power)
Xf = zeros(noTx, Nfft, Nsym);     % [Tx x subcarrier x OFDMsym]
for t = 1:noTx
    if mod == 1
        sym = bpsk_mod_bits(txBits(:,t));
    else
        sym = qpsk_mod_gray_bits(txBits(:,t));
    end
    Xf(t,:,:) = reshape(sym, [Nfft, Nsym]);
end

%% -------------------- OFDM TX (per Tx antenna) --------------------
% IFFT across subcarriers
xt = zeros(noTx, Nfft, Nsym);
for t = 1:noTx
    xt(t,:,:) = ifft(squeeze(Xf(t,:,:)), Nfft, 1);
end

% Add cyclic prefix
xt_cp = zeros(noTx, Nfft+Ncp, Nsym);
for t = 1:noTx
    tmp = squeeze(xt(t,:,:)); % [Nfft x Nsym]
    xt_cp(t,:,:) = [tmp(end-Ncp+1:end,:); tmp];
end

% Serialize each Tx stream
txWave = zeros(noTx, (Nfft+Ncp)*Nsym);
for t = 1:noTx
    txWave(t,:) = reshape(squeeze(xt_cp(t,:,:)), 1, []);
end

%% -------------------- Frequency-Selective 2x2 Rayleigh Channel --------------------
% h_taps: [Lch x noRx x noTx]
h_taps = zeros(Lch, noRx, noTx);
for r = 1:noRx
    for t = 1:noTx
        h_taps(:,r,t) = (randn(Lch,1) + 1i*randn(Lch,1)) .* sqrt(pdp(:)/2);
    end
end

% Frequency response per subcarrier: Hf(:,:,k) is [noRx x noTx]
Hf = zeros(noRx, noTx, Nfft);
for k = 1:Nfft
    % build via FFT of taps for each (r,t)
end
for r = 1:noRx
    for t = 1:noTx
        Hrt = fft([h_taps(:,r,t); zeros(Nfft-Lch,1)], Nfft); % [Nfft x 1]
        Hf(r,t,:) = reshape(Hrt, [1 1 Nfft]);
    end
end

% Pass through channel (time-domain convolution per Rx)
rxClean = zeros(noRx, size(txWave,2) + Lch - 1);
for r = 1:noRx
    y = zeros(1, size(txWave,2) + Lch - 1);
    for t = 1:noTx
        y = y + conv(txWave(t,:), h_taps(:,r,t).');
    end
    rxClean(r,:) = y;
end

% For clean reshaping, take the aligned portion (assume perfect timing)
NsampFrame = (Nfft+Ncp)*Nsym;
rxClean = rxClean(:, 1:NsampFrame); % drop tail (ISI handled by CP if Ncp >= Lch-1)

% Estimate received signal power (for noise scaling)
sigPow = mean(abs(rxClean(:)).^2);

%% -------------------- BER vs SNR Simulation --------------------
ber_zf   = zeros(size(SNR_dB));
ber_sic  = zeros(size(SNR_dB));

for is = 1:numel(SNR_dB)
    % AWGN
    noiseVar = sigPow / SNRlin(is);
    w = sqrt(noiseVar/2) * (randn(size(rxClean)) + 1i*randn(size(rxClean)));
    rxWave = rxClean + w;

    % ---------- OFDM RX ----------
    % Reshape to symbols
    rxCP = reshape(rxWave.', (Nfft+Ncp), Nsym, noRx);    % [time x sym x rx]
    rxCP = permute(rxCP, [3 1 2]);                      % [rx x time x sym]

    % Remove CP
    rxNoCP = rxCP(:, Ncp+1:end, :);                     % [rx x Nfft x sym]

    % FFT to frequency domain
    Yf = zeros(noRx, Nfft, Nsym);
    for r = 1:noRx
        Yf(r,:,:) = fft(squeeze(rxNoCP(r,:,:)), Nfft, 1);
    end

    % ---------- Detect per subcarrier per OFDM symbol ----------
    % Collect detected bits for both receivers
    rxBits_zf  = zeros(Nb, noTx);
    rxBits_sic = zeros(Nb, noTx);

    bitPtr = 1;

    for m = 1:Nsym
        for ksc = 1:Nfft
            yk = squeeze(Yf(:,ksc,m));        % [noRx x 1]
            Hk = squeeze(Hf(:,:,ksc));        % [noRx x noTx]

            % ---- ZF estimate of both streams ----
            % xhat = pinv(Hk)*yk;  % works
            xhat = (Hk' * Hk) \ (Hk' * yk);   % classic ZF (2x2)

            % Decisions for ZF-only
            if mod == 1
                x1_zf = bpsk_decision(xhat(1));
                x2_zf = bpsk_decision(xhat(2));
                bits_zf = [bpsk_demod_bit(x1_zf); bpsk_demod_bit(x2_zf)];
            else
                x1_zf = qpsk_decision_gray(xhat(1));
                x2_zf = qpsk_decision_gray(xhat(2));
                bits_zf = [qpsk_demod_gray_bits(x1_zf); qpsk_demod_gray_bits(x2_zf)];
            end

            % Store bits in stream order: Tx1 then Tx2
            if mod == 1
                rxBits_zf(bitPtr,1) = bits_zf(1);
                rxBits_zf(bitPtr,2) = bits_zf(2);
            else
                rxBits_zf(bitPtr:bitPtr+1,1) = bits_zf(1:2);
                rxBits_zf(bitPtr:bitPtr+1,2) = bits_zf(3:4);
            end

            % ---- ZF-SIC + MRC (match your original logic: detect stream-2 first) ----
            % 1) Detect x2 from ZF
            if mod == 1
                x2_hat = bpsk_decision(xhat(2));
            else
                x2_hat = qpsk_decision_gray(xhat(2));
            end

            % 2) Cancel stream-2: y_res = y - h2*x2_hat
            h2 = Hk(:,2);
            y_res = yk - h2 * x2_hat;

            % 3) MRC detect stream-1 from residual: x1 = (h1^H y_res)/(h1^H h1)
            h1 = Hk(:,1);
            x1_mrc = (h1' * y_res) / (h1' * h1);

            if mod == 1
                x1_hat = bpsk_decision(x1_mrc);
                bits_sic = [bpsk_demod_bit(x1_hat); bpsk_demod_bit(x2_hat)];
                rxBits_sic(bitPtr,1) = bits_sic(1);
                rxBits_sic(bitPtr,2) = bits_sic(2);
            else
                x1_hat = qpsk_decision_gray(x1_mrc);
                bits1 = qpsk_demod_gray_bits(x1_hat);
                bits2 = qpsk_demod_gray_bits(x2_hat);
                rxBits_sic(bitPtr:bitPtr+1,1) = bits1;
                rxBits_sic(bitPtr:bitPtr+1,2) = bits2;
            end

            % advance pointer
            if mod == 1
                bitPtr = bitPtr + 1;
            else
                bitPtr = bitPtr + 2;
            end
        end
    end

    % ---------- BER ----------
    ber_zf(is)  = mean(rxBits_zf(:)  ~= txBits(:));
    ber_sic(is) = mean(rxBits_sic(:) ~= txBits(:));

    fprintf("SNR=%2d dB: BER(ZF)=%.3e, BER(ZF-SIC)=%.3e\n", SNR_dB(is), ber_zf(is), ber_sic(is));
end

%% -------------------- Plot --------------------
figure;
semilogy(SNR_dB, ber_zf,  'o-', 'LineWidth', 2); hold on;
semilogy(SNR_dB, ber_sic, 's-', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
if mod==1
    title(sprintf('2x2 MIMO-OFDM (BPSK), Nfft=%d, Ncp=%d, L=%d, ZF vs ZF-SIC+MRC', Nfft, Ncp, Lch));
else
    title(sprintf('2x2 MIMO-OFDM (QPSK), Nfft=%d, Ncp=%d, L=%d, ZF vs ZF-SIC+MRC', Nfft, Ncp, Lch));
end
legend('ZF-only', 'ZF-SIC + MRC', 'Location','southwest');

%% ===================== Helper Functions =====================

% ---- BPSK ----
function s = bpsk_mod_bits(bits)
    bits = bits(:);
    s = 1 - 2*bits;              % 0->+1, 1->-1
    s = complex(s,0);
end

function x = bpsk_decision(z)
    x = 1;
    if real(z) < 0
        x = -1;
    end
end

function b = bpsk_demod_bit(x)
    % +1 -> 0, -1 -> 1
    b = (real(x) < 0);
end

% ---- QPSK Gray (unit power) ----
function s = qpsk_mod_gray_bits(bits)
    bits = bits(:);
    assert(mod(numel(bits),2)==0, 'QPSK needs even number of bits');
    b = reshape(bits, 2, []).';

    % Gray mapping:
    % 00 -> +1 + j
    % 01 -> -1 + j
    % 11 -> -1 - j
    % 10 -> +1 - j
    s = zeros(size(b,1),1);
    for n = 1:size(b,1)
        if     b(n,1)==0 && b(n,2)==0, s(n)= 1+1i;
        elseif b(n,1)==0 && b(n,2)==1, s(n)=-1+1i;
        elseif b(n,1)==1 && b(n,2)==1, s(n)=-1-1i;
        else                           s(n)= 1-1i; % 10
        end
    end
    s = s / sqrt(2); % unit average power
end

function x = qpsk_decision_gray(z)
    % nearest neighbor to {±1±j}/sqrt(2)
    const = [ 1+1i; -1+1i; -1-1i; 1-1i ] / sqrt(2);
    [~,idx] = min(abs(z-const).^2);
    x = const(idx);
end

function bits = qpsk_demod_gray_bits(x)
    % returns 2 bits (column) for one QPSK symbol (Gray), assumes x is one of constellation points
    x = x * sqrt(2); % undo normalization
    I = real(x); Q = imag(x);
    if I>=0 && Q>=0
        bits = [0;0];
    elseif I<0 && Q>=0
        bits = [0;1];
    elseif I<0 && Q<0
        bits = [1;1];
    else
        bits = [1;0];
    end
end
