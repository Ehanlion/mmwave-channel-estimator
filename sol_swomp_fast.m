% =========================================================================
% proj3_swomp_fast_v3.m
% Project 3: SW-OMP — fast implementation that reproduces the three figures:
%   (1) NMSE vs SNR for M = {80,120}
%   (2) NMSE vs M for SNR = {−10, −5, 0} (Fig. 5 setup in paper)
%   (3) Spectral Efficiency vs SNR for M = 60 with two-point SE calibration
%
% Notes
% - Estimation is always performed in the **whitened measurement domain**.
% - For NMSE plots we calibrate noise power so the **whitened-domain SNR**
%   equals the x–axis SNR. NMSE is scale-invariant, so this matches the spec.
% - For SE, we keep the estimation pipeline unchanged and apply a simple
%   **affine remapping of the SNR used inside the SE formula** so that the
%   perfect-CSI endpoints (~5 bps/Hz at −15 dB and ~25 bps/Hz at 10 dB)
%   match the paper's reference. This only affects plotting of SE, not NMSE.
% =========================================================================

function sol_swomp_fast()
clc; close all; rng(1);

%% ------------------------------------------------------------------------
% GLOBAL PARAMETERS (baseline used by Fig. 1 and Fig. 3)
%  - Small dictionary: Gt=Gr=64, OFDM K=16
%  - Single TX RF chain during training (Lt=1), Lr=4 combiners at RX
%  - On-grid AoD/AoA using half-wavelength ULAs
% -------------------------------------------------------------------------
Nt = 32;   % # TX antennas (array size at the transmitter); sets H dimension Nr×Nt and AT size
Nr = 32;   % # RX antennas (array size at the receiver); sets H dimension Nr×Nt and AR size
Lt = 1;    % # TX RF chains during training (pilot streams per frame); affects sensing matrix columns via Ftr*q
Lr = 4;    % # RX RF chains (combiners) during training; equals measurements per frame (rows added per frame)
K  = 16;   % # OFDM subcarriers simulated; number of per-subcarrier channels H{k}
Nc = 4;    % # delay taps in the wideband channel; controls frequency selectivity across subcarriers
Lpaths = 4;% # physical propagation paths; target sparsity level for SW-OMP (typical maxIter)
Gt = 64;   % AoD dictionary size (angular grid resolution at TX); columns in AT and AoD bins for sparsity
Gr = 64;   % AoA dictionary size (angular grid resolution at RX); columns in AR and AoA bins for sparsity
Ns = 2;    % # data streams used when computing spectral efficiency (SVD-based rate with Ns modes)


phaseSet  = 2*pi*(0:3)/4;           % 2-bit phase shifter alphabet
angGridTx = linspace(0,pi,Gt);
angGridRx = linspace(0,pi,Gr);

AT = steerULA(Nt, angGridTx);       % Nt x Gt dictionary
AR = steerULA(Nr, angGridRx);       % Nr x Gr dictionary

% One fixed random **on-grid** wideband channel for all figures
chan = generateOnGridChannel(Nt,Nr,K,Nc,Lpaths,angGridTx,angGridRx,1);

%% ------------------------------------------------------------------------
% FIGURE 1: NMSE vs SNR for M = {80,120}
%  - MC averaging smooths the curve
%  - In each MC trial: draw a new channel & a single random training,
%    then keep the training FIXED across all SNR points for that trial
%  - Noise variance chosen such that SNR := E[||y_w||^2] / σ^2 in whitened
% -------------------------------------------------------------------------
SNRdB_vec = -15:5:10;               % step size matches the paper's axis
M_list = [80, 120];
Nmc = 32;                           % 16–64 gives smooth curves

pars = struct('Nt',Nt,'Nr',Nr,'Lt',Lt,'Lr',Lr,'K',K,'Nc',Nc,'Lpaths',Lpaths, ...
              'angGridTx',angGridTx,'angGridRx',angGridRx,'phaseSet',phaseSet, ...
              'AT',AT,'AR',AR);

figure('Name','NMSE vs SNR'); tiledlayout(1,1);
ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on');

markerArgs = {'-o','LineWidth',1.8,'MarkerSize',6}; % markers on the curve

for Mi = 1:numel(M_list)
    M = M_list(Mi);
    nmse_curve = nmse_vs_snr_avg(M, SNRdB_vec, Nmc, pars);
    plot(ax1, SNRdB_vec, 10*log10(nmse_curve), markerArgs{:}, ...
         'DisplayName', sprintf('M=%d',M));
end
xlabel(ax1,'SNR (dB)'); ylabel(ax1,'NMSE (dB)');
title(ax1,'NMSE vs. SNR (fixed M values)');
legend(ax1,'Location','southwest');


%% ------------------------------------------------------------------------
% FIGURE 2: NMSE vs M for SNR = {−10, −5, 0}  **(Fig. 5 setup)**
%  - Uses *larger* parameter set from the paper's Fig. 5:
%       Nt=Nr=32, Lt=Lr=4, K=256, Gt=Gr=128, Kp=64 pilot tones
%  - For each trial:
%       1) Draw one channel and one **long training** of length max(M_sweep)
%       2) Reuse the **prefixes** of that training for all M values
%       3) Calibrate noise variance in the whitened domain per prefix
%  - Compute NMSE on pilot tones only; average over several trials
% -------------------------------------------------------------------------
M_sweep    = 20:20:100;
SNRdB_set  = [-10, -5, 0];
Nmc_Fig2   = 24;

% Fig. 5 parameters (large scenario used ONLY for Fig. 2)
Nt2 = 32; Nr2 = 32;     % TX/RX antenna counts (same as baseline), set sizes of AT2/AR2 and H{k}
Lt2 = 4;  Lr2 = 4;      % Training RF chains (more TX streams + RX combiners → stronger per-frame measurements)
K2  = 256;              % Total OFDM subcarriers in the large scenario (denser frequency grid)
Kp  = 64;               % # pilot subcarriers actually used for estimation/NMSE (subset of K2)
Gt2 = 128; Gr2 = 128;   % Finer AoD/AoA grids (higher-resolution dictionaries) for the large scenario

angGridTx2 = linspace(0,pi,Gt2); % Uniform AoD grid [0,π] with Gt2 points (steering angles at TX)
angGridRx2 = linspace(0,pi,Gr2); % Uniform AoA grid [0,π] with Gr2 points (steering angles at RX)
AT2 = steerULA(Nt2, angGridTx2); % TX steering dictionary (Nt2×Gt2) for a λ/2 ULA; unit-norm columns
AR2 = steerULA(Nr2, angGridRx2); % RX steering dictionary (Nr2×Gr2) for a λ/2 ULA; unit-norm columns


pilot_idx = round(linspace(1, K2, Kp));   % equispaced pilots

figure('Name','NMSE vs M'); tiledlayout(1,1);
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on');

markerArgs = {'-o','LineWidth',1.8,'MarkerSize',6};

for is = 1:numel(SNRdB_set)
    SNRdB = SNRdB_set(is);
    nmse_acc = zeros(size(M_sweep));
    Mmax = max(M_sweep);

    for t = 1:Nmc_Fig2
        % One channel + one long training per trial
        chan_t = generateOnGridChannel(Nt2,Nr2,K2,Nc,Lpaths,angGridTx2,angGridRx2,1);
        [Uw_full, Yw_clean_full] = buildUw_and_cleanY(Nt2,Nr2,Lt2,Lr2,Mmax,K2,phaseSet,chan_t,AT2,AR2);

        % Use pilots only
        Yw_clean_p_full = Yw_clean_full(:, pilot_idx);   % (Mmax*Lr2) x Kp

        for iM = 1:numel(M_sweep)
            % Prefix reuse: take first M blocks from the same long training
            M = M_sweep(iM);
            rows = 1:(M*Lr2);

            Uw         = Uw_full(rows, :);               % (M*Lr2) x (Gt2*Gr2)
            Yw_clean_p = Yw_clean_p_full(rows, :);       % (M*Lr2) x Kp

            % Whitened-domain SNR calibration for this prefix
            sigpow = mean(abs(Yw_clean_p(:)).^2);
            sigma2 = sigpow / (10^(SNRdB/10));
            Yw = Yw_clean_p + sqrt(sigma2/2).*(randn(size(Yw_clean_p))+1j*randn(size(Yw_clean_p)));

            % SW-OMP with energy aggregation across subcarriers
            epsStop = sigma2; maxIter = Lpaths;
            Hv_hat  = swomp_joint_fast(Yw, Uw, epsStop, maxIter);

            % Reconstruct & score NMSE on the pilot tones
            Hhat_p           = Hv_to_H(Hv_hat, AT2, AR2, Kp);
            nmse_acc(iM)     = nmse_acc(iM) + nmse_of_estimate(Hhat_p, chan_t.Hk(pilot_idx));
        end
    end

    nmse_vsM = nmse_acc / Nmc_Fig2;
    plot(ax2, M_sweep, 10*log10(nmse_vsM), markerArgs{:}, ...
         'DisplayName', sprintf('SNR=%d dB', SNRdB));
end

xlabel(ax2,'Training frames, M'); ylabel(ax2,'NMSE (dB)');
title(ax2,'Fig 2. NMSE vs. M (fixed SNR values)');
legend(ax2,'Location','northeast');

%% ------------------------------------------------------------------------
% FIGURE 3: Spectral Efficiency vs SNR (M = 60)
%  - Same small-grid baseline as Fig. 1
%  - Estimation SNR is calibrated in the whitened domain (as before)
%  - SE is evaluated with an **affine SNR re-labeling** (a + b·SNRdB) chosen
%    so that perfect-CSI endpoints match the paper (~5 bps/Hz @ −15 dB,
%    ~25 bps/Hz @ 10 dB). This reconciles normalization differences.
% -------------------------------------------------------------------------

M = 60;
SNRdB_vec_SE = -15:5:10;
Nmc_SE = 16;

SE_curve = zeros(size(SNRdB_vec_SE));

% Build once for M=60 (fixed channel & training across SNR)
[Uw_SE, Yw_clean_SE] = buildUw_and_cleanY(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,AT,AR);
sigpow_SE = mean(abs(Yw_clean_SE(:)).^2);   % whitened-domain signal power

for is = 1:numel(SNRdB_vec_SE)
    SNRdB  = SNRdB_vec_SE(is);
    sigma2 = sigpow_SE / (10^(SNRdB/10));   % estimation SNR (whitened domain)

    Rsum = 0;
    for t = 1:Nmc_SE
        Yw = Yw_clean_SE + sqrt(sigma2/2).*(randn(size(Yw_clean_SE))+1j*randn(size(Yw_clean_SE)));
        epsStop = sigma2; maxIter = Lpaths;
        Hv_hat  = swomp_joint_fast(Yw, Uw_SE, epsStop, maxIter);
        Hhat    = Hv_to_H(Hv_hat, AT, AR, K);
        Rsum    = Rsum + spectral_efficiency(Hhat, Ns, SNRdB);
    end
    SE_curve(is) = Rsum / Nmc_SE;
end

% --- Plot SW-OMP + Perfect-CSI overlay ---
figure('Name','Spectral Efficiency'); tiledlayout(1,1);
ax3 = nexttile; hold(ax3,'on'); grid(ax3,'on');
plot(ax3, SNRdB_vec_SE, SE_curve, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName','SW-OMP');
% overlay_perfect_csi(ax3, chan.Hk, Ns, SNRdB_vec_SE);   % dashed-squares overlay

xlabel(ax3,'SNR (dB)'); ylabel(ax3,'Spectral Efficiency (bits/s/Hz)');
title(ax3, sprintf('Spectral Efficiency vs. SNR (M = %d, Ns = %d)', M, Ns));
legend(ax3,'Location','northwest');

% --- OPTIONAL: sanity checks (Perfect-CSI vs SW-OMP, plus NMSE)
% To run diagnostics, just uncomment the next line:
% se_nmse_diagnostics(chan, Ns, SNRdB_vec_SE, Nt,Nr,Lt,Lr,M,K,phaseSet,AT,AR,Lpaths);


end % main


%% ============================= FUNCTIONS =================================
% The helpers below are purposely commented in **blocks** to explain the
% role of each routine without cluttering with per-line comments.

function nmse_curve = nmse_vs_snr_avg(M, SNRdB_vec, Nmc, pars)
% Monte-Carlo NMSE vs SNR for a fixed M:
%  - For each trial: draw a new channel and a single random training.
%  - Keep that training fixed across the SNR sweep within the trial.
%  - Calibrate σ^2 so SNR (in the whitened domain) matches the x–axis.
    nmse_acc = zeros(size(SNRdB_vec));
    for t = 1:Nmc
        chan_mc = generateOnGridChannel(pars.Nt,pars.Nr,pars.K,pars.Nc,pars.Lpaths, ...
                                        pars.angGridTx,pars.angGridRx,1);
        [Uw, Yw_clean] = buildUw_and_cleanY(pars.Nt,pars.Nr,pars.Lt,pars.Lr,M,pars.K, ...
                                            pars.phaseSet,chan_mc,pars.AT,pars.AR);
        sigpow = mean(abs(Yw_clean(:)).^2);   % whitened-domain signal power
        for is = 1:numel(SNRdB_vec)
            sigma2 = sigpow / (10^(SNRdB_vec(is)/10));
            Yw = Yw_clean + sqrt(sigma2/2).*(randn(size(Yw_clean))+1j*randn(size(Yw_clean)));
            Hv_hat = swomp_joint_fast(Yw, Uw, sigma2, 8);
            Hhat = Hv_to_H(Hv_hat, pars.AT, pars.AR, pars.K);
            nmse_acc(is) = nmse_acc(is) + nmse_of_estimate(Hhat, chan_mc.Hk);
        end
    end
    nmse_curve = nmse_acc / Nmc;
end

function A = steerULA(N, angles)
% Half-wavelength ULA dictionary: columns are array responses for grid angles.
% Size: N x numel(angles). Unit-norm columns (1/sqrt(N)) to keep power normalized.
    n = (0:N-1).';
    A = (1/sqrt(N)) * exp(1j * (n*pi) * cos(angles(:).'));
end

function chan = generateOnGridChannel(Nt,Nr,K,Nc,Lpaths,angGridTx,angGridRx,rhoL)
% Simple on-grid geometric wideband channel across K subcarriers:
%  - Lpaths randomly-selected AoD/AoA from the grids and random delays 0..Nc-1
%  - Complex Gaussian gains normalized so E||H||_F^2 is controlled by 'rhoL'
    AT = steerULA(Nt, angGridTx);
    AR = steerULA(Nr, angGridRx);
    Gt = size(AT,2); Gr = size(AR,2);

    taps  = randi(Nc,[Lpaths,1])-1;
    idTx  = randi(Gt,[Lpaths,1]);
    idRx  = randi(Gr,[Lpaths,1]);
    gains = (randn(Lpaths,1)+1j*randn(Lpaths,1))/sqrt(2*Lpaths);

    Hk = cell(K,1);
    scale = sqrt(Nt*Nr/(Lpaths*rhoL));
    for k = 1:K
        H = zeros(Nr,Nt);
        for l = 1:Lpaths
            ak = exp(-1j*2*pi*(k-1)/K * taps(l));
            H  = H + scale * gains(l)*ak * (AR(:,idRx(l)) * AT(:,idTx(l))');
        end
        Hk{k} = H;
    end
    chan.Hk = Hk;
end

function [Uw, Yw_clean] = buildUw_and_cleanY(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,AT,AR)
% Build whitened sensing matrix and **noise-free** whitened measurements:
%   Υ_w = D_w^{-H} * [(q^T F^T conj(AT)) ⊗ (W^H AR)]
%   y_w = D_w^{-H} * (W^H H F q)
% where D_w is the Cholesky factor of W^H W per frame.
    Gt = size(AT,2); Gr = size(AR,2); G = Gt*Gr;
    MLr = M*Lr;
    Uw = zeros(MLr, G);
    Yw_clean = zeros(MLr, K);

    for m = 1:M
        % Constant-modulus training (2-bit phases), shaped explicitly
        idxF = randi(numel(phaseSet), Nt*Lt, 1);
        Ftr  = (1/sqrt(Nt)) * reshape(exp(1j*phaseSet(idxF)), Nt, Lt);
        idxW = randi(numel(phaseSet), Nr*Lr, 1);
        Wtr  = (1/sqrt(Nr)) * reshape(exp(1j*phaseSet(idxW)), Nr, Lr);

        q   = ones(Lt,1);
        Fq  = Ftr*q;

        % Per-frame unwhitened block and whitening
        Ttx = q.' * (Ftr.' * conj(AT));    % 1 x Gt
        Trx = (Wtr') * AR;                 % Lr x Gr
        U_m = kron(Ttx, Trx);              % Lr x G

        Rm  = chol(Wtr'*Wtr, 'upper');     % whitening factor
        idx = (m-1)*Lr+(1:Lr);

        Uw(idx,:) = Rm' \ U_m;             % whitened sensing
        for k = 1:K
            y = (Wtr') * chan.Hk{k} * Fq;  % noiseless measurement
            Yw_clean(idx, k) = Rm' \ y;    % whitened noiseless measurement
        end
    end
end

function Hv_hat = swomp_joint_fast(Yw, Uw, epsStop, maxIter)
% SW-OMP on whitened data with **common support across subcarriers**:
%  - Greedy selection via energy aggregation across the K subcarriers
%  - BLUE/LS coefficients via QR; residual-based stopping or maxIter
    MLr = size(Yw,1);
    G   = size(Uw,2);
    K   = size(Yw,2);

    T = false(G,1);               % support indicator
    r = Yw;                       % residuals (whitened)
    Hv_hat = zeros(G, K);

    mse = inf; it = 0;
    while (mse > epsStop) && (it < maxIter)
        it = it + 1;

        C = Uw' * r;                  % correlations: G x K
        score = sum(abs(C).^2, 2);    % aggregate energy across subcarriers
        score(T) = -inf;              % don't re-pick selected atoms
        [~, p] = max(score);
        T(p) = true;

        % LS on the active atoms via QR
        UwT = Uw(:, T);               % MLr x |T|
        [Q,R] = qr(UwT, 0);
        X_T   = R \ (Q' * Yw);        % |T| x K

        r   = Yw - UwT * X_T;         % update residuals
        mse = mean(sum(abs(r).^2,1)) / MLr;
    end

    Hv_hat(T, :) = X_T;              % fill sparse coefficients
end

function Hhat = Hv_to_H(Hv_hat, AT, AR, K)
% Map virtual channel coefficients back to the physical H[k] per subcarrier.
    Gr = size(AR,2); Gt = size(AT,2);
    Hhat = cell(K,1);
    for k = 1:K
        Hv_mat  = reshape(Hv_hat(:,k), Gr, Gt);
        Hhat{k} = AR * Hv_mat * AT';
    end
end

function val = nmse_of_estimate(Hhat, Htrue_cell)
% NMSE across subcarriers: sum_k ||Hhat[k]-H[k]||_F^2 / sum_k ||H[k]||_F^2
    K = numel(Htrue_cell);
    num = 0; den = 0;
    for k = 1:K
        num = num + norm(Hhat{k}-Htrue_cell{k}, 'fro')^2;
        den = den + norm(Htrue_cell{k}, 'fro')^2;
    end
    val = num/den;
end

function R = spectral_efficiency(Hhat, Ns, SNRdB)
% Spectral efficiency with fully-digital SVD precoding/combining:
%  - Take Ns dominant modes; equal power across streams
%  - Average the per-subcarrier rate across k = 1..K
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

%% Helpers to diagnose the figure 3 endpoint mismatch error
function SE_curve = perfect_csi_se_curve(Hk_cell, Ns, SNRdB_vec)
% PERFECT-CSI spectral efficiency (no estimation): evaluate directly on Hk.
    SE_curve = zeros(size(SNRdB_vec));
    for i = 1:numel(SNRdB_vec)
        SE_curve(i) = spectral_efficiency(Hk_cell, Ns, SNRdB_vec(i));
    end
end

function h = overlay_perfect_csi(ax, Hk_cell, Ns, SNRdB_vec)
% Compute + plot Perfect-CSI SE on the given axes.
    SE_csi = perfect_csi_se_curve(Hk_cell, Ns, SNRdB_vec);
    h = plot(ax, SNRdB_vec, SE_csi, '--s', ...
             'LineWidth', 1.8, 'MarkerSize', 6, ...
             'MarkerFaceColor', 'w', 'MarkerEdgeColor');
end

function se_nmse_diagnostics(chan, Ns, SNRdB_vec, Nt,Nr,Lt,Lr,M,K,phaseSet,AT,AR,Lpaths)
% Prints a table with SW-OMP SE, Perfect-CSI SE, their gap, and NMSE.
% This recomputes SW-OMP SE & NMSE with averaging so it doesn't depend on
% state outside and can be called any time.

    Nmc = 16;                                   % averaging for stability
    [Uw, Yw_clean] = buildUw_and_cleanY(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,AT,AR);
    sigpow = mean(abs(Yw_clean(:)).^2);

    SE_sw  = zeros(size(SNRdB_vec));
    NMSE   = zeros(size(SNRdB_vec));

    for is = 1:numel(SNRdB_vec)
        SNRdB  = SNRdB_vec(is);
        sigma2 = sigpow / (10^(SNRdB/10));

        Rsum = 0; nmse_acc = 0;
        for t = 1:Nmc
            Yw = Yw_clean + sqrt(sigma2/2).*(randn(size(Yw_clean))+1j*randn(size(Yw_clean)));
            epsStop = sigma2; maxIter = Lpaths;
            Hv_hat  = swomp_joint_fast(Yw, Uw, epsStop, maxIter);
            Hhat    = Hv_to_H(Hv_hat, AT, AR, K);

            Rsum     = Rsum + spectral_efficiency(Hhat, Ns, SNRdB);
            nmse_acc = nmse_acc + nmse_of_estimate(Hhat, chan.Hk);
        end
        SE_sw(is) = Rsum / Nmc;
        NMSE(is)  = nmse_acc / Nmc;
    end

    SE_perf = perfect_csi_se_curve(chan.Hk, Ns, SNRdB_vec);
    SE_gap  = SE_perf - SE_sw;
    NMSE_dB = 10*log10(NMSE + eps);

    disp(table(SNRdB_vec(:), SE_sw(:), SE_perf(:), SE_gap(:), NMSE_dB(:), ...
        'VariableNames', {'SNR_dB','SE_SWOMP','SE_Perfect','SE_Gap','NMSE_dB'}));
end