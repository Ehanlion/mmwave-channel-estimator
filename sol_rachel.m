%% ECE 233 – Project 3 (Compressive mmWave Channel Estimation, on-grid)
% Reproduces the three figures requested in the project (NMSE vs SNR for M=80/120,
% NMSE vs M at SNR={-10,-5,0} dB, and Spectral Efficiency vs SNR for M=60).
% Reference model and metrics follow [R-4] (Rodríguez-Fernández et al., TWC 2018).

clear; clc; rng(1);

%% ------------------------- 1) Global parameters --------------------------
Nt = 32; Nr = 32;              % antennas (Tx,Rx)
Lt = 1;  Lr = 4;               % RF chains (training)
K  = 16;                       % OFDM subcarriers (paper's baseline)
Nc = 4;                        % delay taps
Lpaths = 4;                    % number of physical paths
Gt = 64; Gr = 64;              % dictionary sizes (AoD, AoA grid)
Ns = 2;                        % streams used for spectral efficiency (post-estimation)
rhoL = 1;                      % path loss normalization
Ptx = 1;                       % total TX power

phaseSet = 2*pi*(0:3)/4;       % 2-bit phases for training beams
angGridTx = linspace(0,pi,Gt); % on-grid AoD grid
angGridRx = linspace(0,pi,Gr); % on-grid AoA grid

% Steering-vector dictionaries (half-wavelength ULA)
AT = steerULA(Nt, angGridTx);
AR = steerULA(Nr, angGridRx);

% Kronecker dictionary Ψ = (conj(AT) ⊗ AR)  [size (Nt*Nr) x (Gt*Gr)]
Psi = kron(conj(AT), AR);

% Pre-generate on-grid wideband channel (fixed angles chosen from the grids)
chan = generateOnGridChannel(Nt,Nr,K,Nc,Lpaths,angGridTx,angGridRx,rhoL);

%% ------------------------- 2) Figure 1: NMSE vs SNR ----------------------
SNRdB_vec = -15:2.5:10;                 % SNR range (dB)
M_list = [80, 120];                     % training frames

figure('Name','NMSE vs SNR'); tiledlayout(1,1); ax1=nexttile; hold on; grid on;
leg = {};

for Mi = 1:numel(M_list)
    M = M_list(Mi);
    nmse_curve = zeros(size(SNRdB_vec));
    for is = 1:numel(SNRdB_vec)
        SNRdB = SNRdB_vec(is);
        sigma2 = Ptx/(10^(SNRdB/10));  % SNR = P / sigma^2 (ρ_L=1) [R-4]
        % ---------------- Training (random per-frame beams) ----------------
        [Phi, Y_all, Dwinv] = buildMeasurements(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,sigma2);
        % ---------------- SW-OMP (whitened, common support) ---------------
        epsStop = sigma2; maxIter = 8;  % up to 8 paths; stop when MSE <= σ^2
        Hv_hat = swomp_joint(Y_all, Phi, Psi, Dwinv, K, epsStop, maxIter);
        % ---------------- Reconstruct H[k] and compute NMSE ---------------
        Hhat = reshapeHvToH(Hv_hat, AT, AR, K);
        nmse_curve(is) = nmse_of_estimate(Hhat, chan.Hk);
    end
    plot(ax1, SNRdB_vec, 10*log10(nmse_curve), 'LineWidth', 1.8, 'DisplayName',sprintf('M=%d',M));
    leg{end+1} = sprintf('M=%d',M);
end
xlabel('SNR (dB)'); ylabel('NMSE (dB)'); title('NMSE vs SNR (on-grid AoA/AoD)');
legend(leg, 'Location','southwest');

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
        [Phi, Y_all, Dwinv] = buildMeasurements(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,sigma2);
        epsStop = sigma2; maxIter = 8;
        Hv_hat = swomp_joint(Y_all, Phi, Psi, Dwinv, K, epsStop, maxIter);
        Hhat = reshapeHvToH(Hv_hat, AT, AR, K);
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

for is = 1:numel(SNRdB_vec_SE)
    SNRdB = SNRdB_vec_SE(is);
    sigma2 = Ptx/(10^(SNRdB/10));
    [Phi, Y_all, Dwinv] = buildMeasurements(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,sigma2);
    epsStop = sigma2; maxIter = 8;
    Hv_hat = swomp_joint(Y_all, Phi, Psi, Dwinv, K, epsStop, maxIter);
    Hhat = reshapeHvToH(Hv_hat, AT, AR, K);
    SE_curve(is) = spectral_efficiency(Hhat, Ns, SNRdB);
end

figure('Name','Spectral Efficiency'); plot(SNRdB_vec_SE, SE_curve,'-o','LineWidth',1.8);
grid on; xlabel('SNR (dB)'); ylabel('Spectral Efficiency (bits/s/Hz)');
title(sprintf('Spectral Efficiency vs SNR (M = %d, Ns = %d)', M, Ns));

%% ============================= FUNCTIONS ================================

function A = steerULA(N, angles)
% Half-wavelength ULA steering dictionary, size N x numel(angles)
    n = (0:N-1).';
    A = (1/sqrt(N)) * exp(1j * (n*pi) .* cos(angles(:).'));
end

function chan = generateOnGridChannel(Nt,Nr,K,Nc,Lpaths,angGridTx,angGridRx,rhoL)
% On-grid wideband geometric channel per [R-4, Sec. II].
% Randomly select Lpaths indices on the AoD/AoA grids; random delays 0..Nc-1
    Gt = numel(angGridTx); Gr = numel(angGridRx);
    AT = steerULA(Nt, angGridTx);
    AR = steerULA(Nr, angGridRx);

    taps = randi(Nc,[Lpaths,1])-1;
    idTx = randi(Gt,[Lpaths,1]);
    idRx = randi(Gr,[Lpaths,1]);
    gains = (randn(Lpaths,1)+1j*randn(Lpaths,1))/sqrt(2*Lpaths); % normalize

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

function [Phi, Y_all, Dwinv] = buildMeasurements(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,sigma2)
% Build stacked sensing matrix Φ, whiteners, and stacked measurements y[k]
% y[k] has length MLr, Φ is (MLr) x (NtNr).  Whitening uses Dw = blkdiag(Wm'*Wm)^(1/2).
    NtNr = Nt*Nr;
    Y_all = zeros(M*Lr, K);    % columns are y[k]
    Phi_blocks = cell(M,1);
    Dblocks = cell(M,1);

    for m = 1:M
        % Random fully-connected training precoder/combiner (quantized phases)
        Ftr = (1/sqrt(Nt)) * exp(1j*phaseSet(randi(numel(phaseSet), Nt, Lt)));
        Wtr = (1/sqrt(Nr)) * exp(1j*phaseSet(randi(numel(phaseSet), Nr, Lr)));
        q   = ones(Lt,1);                      % frequency-flat Lt×1
        % Per-frame measurement matrix Φ(m) = (q^T Ftr^T ⊗ Wtr^*)
        Phi_m = kron( (q.'*Ftr.').', Wtr' );    % size Lr x (NtNr)
        Phi_blocks{m} = Phi_m;

        % Collect measurements across subcarriers
        for k = 1:K
            ymk = (Wtr') * chan.Hk{k} * (Ftr*q) + sqrt(sigma2/2) * (randn(Lr,1)+1j*randn(Lr,1));
            Y_all((m-1)*Lr+(1:Lr), k) = ymk;
        end

        % Whitening block for this frame: Dw_block = chol(Wtr'*Wtr)
        Dw_block = chol(Wtr'*Wtr, 'upper');   % Lr x Lr, upper-tri
        Dblocks{m} = Dw_block;
    end

    Phi = vertcat(Phi_blocks{:});                         % (M*Lr) x (NtNr)

    % Assemble block-diagonal Dw and its inverse (for whitening)
    Dw = blkdiag(Dblocks{:});                             % (M*Lr) x (M*Lr)
    Dwinv = inv(Dw);                                      % use small Lr so direct inv OK

    % Apply whitening to measurements once here (y_w = Dwinv' * y)
    Y_all = (Dwinv.') * Y_all;                            % whitened y[k] columns
    % (We pass Dwinv and Φ; the SW-OMP function will whiten Φ the same way.)
end

function Hv_hat = swomp_joint(Y_all, Phi, Psi, Dwinv, K, epsStop, maxIter)
% Simple SW-OMP (whitened) with common support across subcarriers.
% Inputs:
%   Y_all : (MLr x K) matrix; each column is y[k] already whitened (Dwinv' applied)
%   Phi   : (MLr x NtNr) unwhitened sensing matrix
%   Psi   : (NtNr x GtGr) dictionary
%   Dwinv : (MLr x MLr) inverse Cholesky from whitening
% Output:
%   Hv_hat : (GtGr x K) estimated sparse virtual channel per subcarrier.

    MLr = size(Phi,1);
    G = size(Psi,2);
    % Whiten sensing matrix once: Υ_w = Dwinv' * Φ * Ψ  (matches Alg.1 lines 2–3)
    Uw = (Dwinv.') * (Phi * Psi);                         % (MLr x G)

    T = [];                                               % support set (indices in 1..G)
    r = Y_all;                                            % residuals (already whitened)
    Hv_hat = zeros(G, K);                                 % sparse coeffs

    mse = inf; it = 0;
    while (mse > epsStop) && (it < maxIter)
        it = it + 1;
        % Correlations for all k: c[k] = Υ_w^H * r[k]
        C = Uw' * r;                                      % (G x K)
        % Aggregate across subcarriers to exploit common support
        score = sum(abs(C), 2);                           % (G x 1)
        % Choose best new atom not already selected
        [~, p] = max(score .* (~ismember(1:G, T)).');     % scalar index
        T = [T, p];                                       %#ok<AGROW>

        % Weighted LS (BLUE): x_T[k] = (Υ_w(:,T))^\ y_w[k], for all k
        UwT = Uw(:, T);                                   % (MLr x |T|)
        X_T = UwT \ Y_all;                                % (|T| x K), least-squares
        % Update residuals
        r = Y_all - UwT * X_T;                            % (MLr x K)
        mse = mean(sum(abs(r).^2,1)) / MLr;               % average per-sample energy
    end

    % Fill sparse coeffs
    Hv_hat(T, :) = X_T;
end

function Hhat = reshapeHvToH(Hv_hat, AT, AR, K)
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
% NMSE across subcarriers (eq. (34) in [R-4])
    K = numel(Htrue_cell);
    num = 0; den = 0;
    for k = 1:K
        num = num + norm(Hhat{k}-Htrue_cell{k}, 'fro')^2;
        den = den + norm(Htrue_cell{k}, 'fro')^2;
    end
    val = num/den;
end

function R = spectral_efficiency(Hhat, Ns, SNRdB)
% Spectral efficiency with fully-digital SVD TX/RX (eq. (36) in [R-4])
    K = numel(Hhat);
    SNRlin = 10^(SNRdB/10);
    R = 0;
    for k = 1:K
        [U,S,V] = svd(Hhat{k});
        Ue = U(:,1:Ns); Ve = V(:,1:Ns);
        Heff = Ue' * Hhat{k} * Ve;
        s = svd(Heff);
        R = R + sum(log2(1 + (SNRlin/Ns) * (s.^2)));
    end
    R = R / K;
end
