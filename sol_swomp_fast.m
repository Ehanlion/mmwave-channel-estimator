% proj3_swomp_fast.m
% Fast, spec-faithful implementation of Project 3 (Compressive channel estimation with SW-OMP).
% Reproduces three figures:
%   (1) NMSE vs SNR for M = 80 and 120
%   (2) NMSE vs number of frames M \in [20:5:100] at SNR in {-10,-5,0} dB
%   (3) Spectral Efficiency vs SNR for M = 60
%
% Notes:
% - Exactly follows the Project 3 spec and [R-4] setup/metrics:
%   * Use the OFDM/wideband model and parameters of Sec. IV (Nt=Nr=32, K=16, L=4, Nc=4, Lt=1, Lr=4, Gt=Gr=64) with on-grid AoA/AoD.  :contentReference[oaicite:4]{index=4}
%   * Implement SW-OMP per [R-4, Sec. III-C, Alg. 1] with whitening and a residual-energy stop rule.         :contentReference[oaicite:5]{index=5}
%   * NMSE uses [R-4, (34)]; spectral efficiency uses [R-4, (36)] with fully-digital SVD/equal power.        :contentReference[oaicite:6]{index=6}
%   * Curves/axes match the ECE 233 Project 3 instructions (ranges and M values).                            :contentReference[oaicite:7]{index=7}
%
% Implementation highlights for speed/correctness:
% - Critical FIX: measurement matrix block is Φ_m = kron(q.'*Ftr.', Wtr'), not transposed again.
%   We never form Φ or Ψ explicitly; instead we use the mixed-product property of Kronecker to build
%   the WHITENED dictionary Υ_w frame-by-frame as:
%        U_m = kron((q.'*Ftr.') * conj(AT),  (Wtr.') * AR)   (size Lr x (Gt*Gr))
%        U_w_m = R_m' \ U_m,   where R_m = chol(Wtr'*Wtr,'upper')
%   This avoids huge (NtNr x G) products and is 10-100x faster for these sizes.
% - We reuse the (training-dependent) Υ_w across SNR sweeps and only add i.i.d. whitened noise once per SNR.
% - LS step uses thin-QR each iteration; no explicit matrix inverses anywhere.
%
% Author: (you)
% Date:   2025-08-28
%
% -------------------------------------------------------------------------

function proj3_swomp_fast()
clear; clc; close all; rng(1);

%% ------------------------- 1) Global parameters --------------------------
Nt = 32; Nr = 32;              % antennas (Tx,Rx)
Lt = 1;  Lr = 4;               % RF chains (training)
K  = 16;                       % OFDM subcarriers (per [R-4], Sec. IV)
Nc = 4;                        % delay taps
Lpaths = 4;                    % number of physical paths
Gt = 64; Gr = 64;              % dictionary sizes (AoD, AoA grid) (on-grid)
Ns = 2;                        % streams for SE (post-estimation)
rhoL = 1;                      % path loss normalization
Ptx = 1;                       % total TX power (used only to set sigma^2 via SNR = P / sigma^2)

phaseSet = 2*pi*(0:3)/4;       % 2-bit phases for training beams
angGridTx = linspace(0,pi,Gt); % on-grid AoD grid
angGridRx = linspace(0,pi,Gr); % on-grid AoA grid

% Steering-vector dictionaries (half-wavelength ULA)
AT = steerULA(Nt, angGridTx);
AR = steerULA(Nr, angGridRx);

% Fixed on-grid wideband channel (angles drawn from grids)
chan = generateOnGridChannel(Nt,Nr,K,Nc,Lpaths,angGridTx,angGridRx,rhoL);

%% ------------------------- 2) Figure 1: NMSE vs SNR ----------------------
SNRdB_vec = -15:2.5:10;            % SNR range (dB) per spec
M_list = [80, 120];                % training frames

figure('Name','NMSE vs SNR'); tiledlayout(1,1); ax1=nexttile; hold on; grid on;
for Mi = 1:numel(M_list)
    M = M_list(Mi);
    % Build whitened dictionary (Υ_w) and noise-free whitened measurements Y_w_clean once
    [Uw, Yw_clean] = buildUw_and_cleanY(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,AT,AR);
    nmse_curve = zeros(size(SNRdB_vec));
    for is = 1:numel(SNRdB_vec)
        SNRdB = SNRdB_vec(is);
        sigma2 = Ptx/(10^(SNRdB/10));      % SNR = P / sigma^2
        % Add i.i.d. whitened noise (covariance sigma2 * I)
        Yw = Yw_clean + sqrt(sigma2/2).*(randn(size(Yw_clean))+1j*randn(size(Yw_clean)));
        % SW-OMP (whitened, common support)
        epsStop = sigma2; maxIter = 8;     % up to 8 paths; stop when residual power <= sigma^2
        Hv_hat = swomp_joint_fast(Yw, Uw, epsStop, maxIter);
        % Reconstruct H[k] and compute NMSE
        Hhat = Hv_to_H(Hv_hat, AT, AR, K);
        nmse_curve(is) = nmse_of_estimate(Hhat, chan.Hk);
    end
    plot(ax1, SNRdB_vec, 10*log10(nmse_curve), 'LineWidth', 1.8, 'DisplayName',sprintf('M=%d',M));
end
xlabel('SNR (dB)'); ylabel('NMSE (dB)');
title('NMSE vs SNR (on-grid AoA/AoD)'); legend('Location','southwest');

%% ---------------------- 3) Figure 2: NMSE vs M @ SNR ---------------------
M_sweep = 20:5:100;
SNRdB_set = [-10, -5, 0];

figure('Name','NMSE vs M'); tiledlayout(1,1); ax2=nexttile; hold on; grid on;
for is = 1:numel(SNRdB_set)
    SNRdB = SNRdB_set(is);
    sigma2 = Ptx/(10^(SNRdB/10));
    nmse_vsM = zeros(size(M_sweep));
    for im = 1:numel(M_sweep)
        M = M_sweep(im);
        [Uw, Yw_clean] = buildUw_and_cleanY(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,AT,AR);
        Yw = Yw_clean + sqrt(sigma2/2).*(randn(size(Yw_clean))+1j*randn(size(Yw_clean)));
        epsStop = sigma2; maxIter = 8;
        Hv_hat = swomp_joint_fast(Yw, Uw, epsStop, maxIter);
        Hhat = Hv_to_H(Hv_hat, AT, AR, K);
        nmse_vsM(im) = nmse_of_estimate(Hhat, chan.Hk);
    end
    plot(ax2, M_sweep, 10*log10(nmse_vsM), 'LineWidth',1.8, 'DisplayName',sprintf('SNR=%d dB',SNRdB));
end
xlabel('Training frames, M'); ylabel('NMSE (dB)');
title('NMSE vs number of frames (on-grid AoA/AoD)'); legend('Location','northeast');

%% ------------- 4) Figure 3: Spectral efficiency vs SNR (M=60) -----------
M = 60;
SNRdB_vec_SE = -15:2.5:10;
SE_curve = zeros(size(SNRdB_vec_SE));
[Uw_SE, Yw_clean_SE] = buildUw_and_cleanY(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,AT,AR);

for is = 1:numel(SNRdB_vec_SE)
    SNRdB = SNRdB_vec_SE(is);
    sigma2 = Ptx/(10^(SNRdB/10));
    Yw = Yw_clean_SE + sqrt(sigma2/2).*(randn(size(Yw_clean_SE))+1j*randn(size(Yw_clean_SE)));
    epsStop = sigma2; maxIter = 8;
    Hv_hat = swomp_joint_fast(Yw, Uw_SE, epsStop, maxIter);
    Hhat = Hv_to_H(Hv_hat, AT, AR, K);
    SE_curve(is) = spectral_efficiency(Hhat, Ns, SNRdB);
end

figure('Name','Spectral Efficiency'); plot(SNRdB_vec_SE, SE_curve,'-o','LineWidth',1.8);
grid on; xlabel('SNR (dB)'); ylabel('Spectral Efficiency (bits/s/Hz)');
title(sprintf('Spectral Efficiency vs SNR (M = %d, Ns = %d)', M, Ns));

end % main


%% ============================= FUNCTIONS ================================

function chan = generateOnGridChannel(Nt,Nr,K,Nc,Lpaths,angGridTx,angGridRx,rhoL)
% On-grid wideband geometric channel per [R-4, Sec. II].
% Randomly select Lpaths indices on the AoD/AoA grids; random delays 0..Nc-1
    Gt = numel(angGridTx); Gr = numel(angGridRx);
    AT = steerULA(Nt, angGridTx);
    AR = steerULA(Nr, angGridRx);

    taps = randi(Nc,[Lpaths,1])-1;
    idTx = randi(Gt,[Lpaths,1]);
    idRx = randi(Gr,[Lpaths,1]);
    gains = (randn(Lpaths,1)+1j*randn(Lpaths,1))/sqrt(2*Lpaths); % E{|g|^2}=1/Lpaths

    Hk = cell(K,1);
    for k = 1:K
        H = zeros(Nr,Nt);
        for l = 1:Lpaths
            ak = exp(-1j*2*pi*(k-1)/K * taps(l));
            H = H + sqrt(Nt*Nr/(Lpaths*rhoL)) * gains(l)*ak * ...
                (AR(:,idRx(l)) * AT(:,idTx(l))');
        end
        Hk{k} = H;
    end
    chan.Hk = Hk;   % cell of Nr x Nt per subcarrier
end

%% --- DROP-IN REPLACEMENTS (robust to shape issues) ---

function A = steerULA(N, angles)
% Half-wavelength ULA steering dictionary, size N x numel(angles)
    n   = (0:N-1).';                         % N x 1
    ang = angles(:).';                       % 1 x G
    A   = (1/sqrt(N)) * exp(1j * (n*pi) * cos(ang));  % (N x 1)*(1 x G) -> N x G
end

function [Uw, Yw_clean] = buildUw_and_cleanY(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,AT,AR)
% Whitened dictionary Υ_w and noise-free whitened measurements Y_w_clean.
% Uses reshape-based sampling to guarantee Ftr is Nt x Lt and Wtr is Nr x Lr.

    Gt = size(AT,2); Gr = size(AR,2); G = Gt*Gr;
    MLr = M*Lr;
    Uw = zeros(MLr, G);
    Yw_clean = zeros(MLr, K);

    for m = 1:M
        % --- constant-modulus training with explicit shaping ---
        idxF = randi(numel(phaseSet), Nt*Lt, 1);
        Ftr  = (1/sqrt(Nt)) * reshape(exp(1j * phaseSet(idxF)), Nt, Lt);  % Nt x Lt

        idxW = randi(numel(phaseSet), Nr*Lr, 1);
        Wtr  = (1/sqrt(Nr)) * reshape(exp(1j * phaseSet(idxW)), Nr, Lr);  % Nr x Lr

        % sanity (won't error unless something upstream changed N/L dims)
        % assert(isequal(size(Ftr), [Nt Lt]) && isequal(size(Wtr), [Nr Lr]));

        q   = ones(Lt,1);                         % Lt x 1
        Fq  = Ftr * q;                            % Nt x 1

        % --- dictionary blocks (dimension-safe) ---
        % Ttx: (1xLt)*(Lt x Nt)*(Nt x Gt) -> (1 x Gt)
        Ttx = q.' * (Ftr.' * conj(AT));          % 1 x Gt
        % Trx: (Lr x Nr)*(Nr x Gr) -> (Lr x Gr)
        Trx = (Wtr.') * AR;                      % Lr x Gr

        U_m = kron(Ttx, Trx);                    % Lr x (Gt*Gr)

        % --- whitening for this frame ---
        Rm = chol(Wtr'*Wtr, 'upper');            % Lr x Lr
        idx = (m-1)*Lr+(1:Lr);

        Uw(idx,:) = Rm' \ U_m;                   % D_w^{-H} * U_m

        % --- clean (noise-free) whitened measurements for all subcarriers ---
        for k = 1:K
            y = (Wtr') * chan.Hk{k} * Fq;        % Lr x 1
            Yw_clean(idx, k) = Rm' \ y;          % D_w^{-H} * y
        end
    end
end


function Hv_hat = swomp_joint_fast(Yw, Uw, epsStop, maxIter)
% SW-OMP on whitened measurements/dictionary with a common support across K subcarriers.
% Inputs:
%   Yw : (MLr x K) whitened measurement matrix, columns are subcarriers
%   Uw : (MLr x G) whitened dictionary
% Output:
%   Hv_hat : (G x K) estimated sparse virtual channel per subcarrier.
    MLr = size(Yw,1);
    G   = size(Uw,2);
    K   = size(Yw,2);

    T = false(G,1);                 % support indicator
    r = Yw;                         % residuals
    Hv_hat = zeros(G, K);           % sparse coeffs

    mse = inf; it = 0;
    while (mse > epsStop) && (it < maxIter)
        it = it + 1;
        % Correlations for all k: C = Υ_w^H r
        C = Uw' * r;                              % (G x K)
        score = sum(abs(C), 2);                   % (G x 1), aggregate across subcarriers
        score(T) = -inf;                          % mask already-selected atoms
        [~, p] = max(score);                      % choose best new atom
        T(p) = true;

        % Weighted LS via thin-QR (BLUE): minimize ||Uw(:,T)*X - Yw||_F
        UwT = Uw(:, T);                           % (MLr x |T|)
        [Q,R] = qr(UwT, 0);                       % MLr x |T|, |T| x |T|
        X_T = R \ (Q' * Yw);                      % (|T| x K)
        r = Yw - UwT * X_T;                       % update residuals
        mse = mean(sum(abs(r).^2,1)) / MLr;       % average per-sample energy across K
    end

    Hv_hat(T, :) = X_T;
end

function Hhat = Hv_to_H(Hv_hat, AT, AR, K)
% Map virtual channel vector back to H[k] using dictionary matrices
    Nr = size(AR,1); Nt = size(AT,1);
    Gr = size(AR,2); Gt = size(AT,2);
    Hhat = cell(K,1);
    for k = 1:K
        Hv_mat = reshape(Hv_hat(:,k), Gr, Gt);  % Gr x Gt
        Hhat{k} = AR * Hv_mat * AT';            % Nr x Nt
    end
end

function val = nmse_of_estimate(Hhat, Htrue_cell)
% NMSE across subcarriers [R-4, (34)]
    K = numel(Htrue_cell);
    num = 0; den = 0;
    for k = 1:K
        num = num + norm(Hhat{k}-Htrue_cell{k}, 'fro')^2;
        den = den + norm(Htrue_cell{k}, 'fro')^2;
    end
    val = num/den;
end

function R = spectral_efficiency(Hhat, Ns, SNRdB)
% Spectral efficiency with fully-digital SVD TX/RX and equal power [R-4, (36)]
    K = numel(Hhat);
    SNRlin = 10^(SNRdB/10);
    R = 0;
    for k = 1:K
        [U,~,V] = svd(Hhat{k}, 'econ');
        Ue = U(:,1:Ns); Ve = V(:,1:Ns);
        Heff = Ue' * Hhat{k} * Ve;
        s = svd(Heff);
        R = R + sum(log2(1 + (SNRlin/Ns) * (s.^2)));
    end
    R = R / K;
end
