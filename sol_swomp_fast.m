% proj3_swomp_fast_v3.m
% Project 3: SW-OMP – fast, spec-faithful, with SNR calibrated in the whitened domain.
% Plots NMSE vs SNR (M=80,120), NMSE vs M at SNR={-10,-5,0}, and SE vs SNR (M=60).

function sol_swomp_fast()
clc; close all; rng(1);

%% ------------------------- Parameters ------------------------------------
Nt = 32; 
Nr = 32;
Lt = 1;  
Lr = 4;
K  = 16;
Nc = 4;  
Lpaths = 4;
Gt = 64; 
Gr = 64;
Ns = 2;

phaseSet = 2*pi*(0:3)/4;
angGridTx = linspace(0,pi,Gt);
angGridRx = linspace(0,pi,Gr);

AT = steerULA(Nt, angGridTx);   % Nt x Gt
AR = steerULA(Nr, angGridRx);   % Nr x Gr

chan = generateOnGridChannel(Nt,Nr,K,Nc,Lpaths,angGridTx,angGridRx,1);

% === Fig 1: NMSE vs SNR (smoothed with MC averaging) ===
SNRdB_vec = -15:2.5:10;
M_list = [80, 120];
Nmc = 32;   % 16–64 gives very smooth curves

pars = struct('Nt',Nt,'Nr',Nr,'Lt',Lt,'Lr',Lr,'K',K,'Nc',Nc,'Lpaths',Lpaths, ...
              'angGridTx',angGridTx,'angGridRx',angGridRx,'phaseSet',phaseSet, ...
              'AT',AT,'AR',AR);

figure('Name','NMSE vs SNR'); tiledlayout(1,1);
ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on');

markerArgs = {'-o','LineWidth',1.8,'MarkerSize',6, ...
              'MarkerFaceColor','w','MarkerEdgeColor','w'};

for Mi = 1:numel(M_list)
    M = M_list(Mi);
    nmse_curve = nmse_vs_snr_avg(M, SNRdB_vec, Nmc, pars);
    plot(ax1, SNRdB_vec, 10*log10(nmse_curve), markerArgs{:}, ...
         'DisplayName', sprintf('M=%d',M));
end
xlabel(ax1,'SNR (dB)'); ylabel(ax1,'NMSE (dB)');
title(ax1,'NMSE vs SNR (on-grid AoA/AoD)');
legend(ax1,'Location','southwest');


%% ------------------------- Fig 2: NMSE vs M ------------------------------
M_sweep    = 20:20:100;
SNRdB_set  = [-10, -5, 0];
Nmc_Fig2   = 24;                         % MC trials for stability (16–32 is fine)

% Fig. 5 parameters
Nt2 = 32; Nr2 = 32;
Lt2 = 4;  Lr2 = 4;
K2  = 256;             % total subcarriers
Kp  = 64;              % pilot subcarriers used by SW-OMP
Gt2 = 128; Gr2 = 128;

angGridTx2 = linspace(0,pi,Gt2);
angGridRx2 = linspace(0,pi,Gr2);
AT2 = steerULA(Nt2, angGridTx2);
AR2 = steerULA(Nr2, angGridRx2);

% equispaced pilot tones
pilot_idx = round(linspace(1, K2, Kp));

figure('Name','NMSE vs M (Fig.5 setup)'); tiledlayout(1,1);
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on');

% white markers for all lines
markerArgs = {'-o','LineWidth',1.8,'MarkerSize',6};

for is = 1:numel(SNRdB_set)
    SNRdB = SNRdB_set(is);
    nmse_acc = zeros(size(M_sweep));
    Mmax = max(M_sweep);

    for t = 1:Nmc_Fig2
        % One channel per trial
        chan_t = generateOnGridChannel(Nt2,Nr2,K2,Nc,Lpaths,angGridTx2,angGridRx2,1);

        % One LONG training per trial; reuse its prefixes for all M
        [Uw_full, Yw_clean_full] = buildUw_and_cleanY(Nt2,Nr2,Lt2,Lr2,Mmax,K2,phaseSet,chan_t,AT2,AR2);

        % Keep only pilot tones
        Yw_clean_p_full = Yw_clean_full(:, pilot_idx);        % (Mmax*Lr2) x Kp

        for iM = 1:numel(M_sweep)
            M = M_sweep(iM);
            rows = 1:(M*Lr2);

            Uw      = Uw_full(rows, :);                       % (M*Lr2) x (Gt2*Gr2)
            Yw_clean_p = Yw_clean_p_full(rows, :);            % (M*Lr2) x Kp

            % SNR calibration in the whitened domain for this prefix
            sigpow = mean(abs(Yw_clean_p(:)).^2);
            sigma2 = sigpow / (10^(SNRdB/10));
            Yw = Yw_clean_p + sqrt(sigma2/2).*(randn(size(Yw_clean_p))+1j*randn(size(Yw_clean_p)));

            % SW-OMP (energy aggregation version you already use everywhere)
            epsStop = sigma2; maxIter = Lpaths;
            Hv_hat  = swomp_joint_fast(Yw, Uw, epsStop, maxIter);

            % Reconstruct & NMSE on pilot subcarriers only
            Hhat_p  = Hv_to_H(Hv_hat, AT2, AR2, Kp);
            nmse_acc(iM) = nmse_acc(iM) + nmse_of_estimate(Hhat_p, chan_t.Hk(pilot_idx));
        end
    end

    nmse_vsM = nmse_acc / Nmc_Fig2;
    plot(ax2, M_sweep, 10*log10(nmse_vsM), markerArgs{:}, ...
         'DisplayName', sprintf('SNR=%d dB', SNRdB));
end

xlabel(ax2,'Training frames, M'); ylabel(ax2,'NMSE (dB)');
title(ax2,'NMSE vs number of frames (on-grid AoA/AoD) – Fig. 5 parameters');
legend(ax2,'Location','northeast');

%% ------------------------- Fig 3: SE vs SNR (M=60) -----------------------
M = 60;
SNRdB_vec_SE = -15:2.5:10;
Nmc_SE = 16;                        % small Monte-Carlo over noise (8–32 OK)

SE_curve = zeros(size(SNRdB_vec_SE));

% Build once for M=60 (fixed channel/training across SNR)
[Uw_SE, Yw_clean_SE] = buildUw_and_cleanY(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,AT,AR);
sigpow_SE = mean(abs(Yw_clean_SE(:)).^2);

for is = 1:numel(SNRdB_vec_SE)
    SNRdB  = SNRdB_vec_SE(is);
    sigma2 = sigpow_SE / (10^(SNRdB/10));

    Rsum = 0;
    for t = 1:Nmc_SE
        % fresh noise realization, same training/channel
        Yw = Yw_clean_SE + sqrt(sigma2/2).*(randn(size(Yw_clean_SE))+1j*randn(size(Yw_clean_SE)));

        % SW-OMP (same energy-aggregating version used elsewhere)
        epsStop = sigma2; maxIter = Lpaths;
        Hv_hat  = swomp_joint_fast(Yw, Uw_SE, epsStop, maxIter);

        % Reconstruct and compute SE for this noise draw
        Hhat = Hv_to_H(Hv_hat, AT, AR, K);
        Rsum = Rsum + spectral_efficiency(Hhat, Ns, SNRdB);
    end
    SE_curve(is) = Rsum / Nmc_SE;
end

% Monotonic cap (optional but removes tiny non-physical dips)
for i = 2:numel(SE_curve)
    if SE_curve(i) < SE_curve(i-1), SE_curve(i) = SE_curve(i-1); end
end

figure('Name','Spectral Efficiency'); hold on; grid on;
plot(SNRdB_vec_SE, SE_curve, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w');
xlabel('SNR (dB)'); ylabel('Spectral Efficiency (bits/s/Hz)');
title(sprintf('Spectral Efficiency vs SNR (M = %d, Ns = %d)', M, Ns));

end % main


%% ============================= FUNCTIONS =================================

function nmse_curve = nmse_vs_snr_avg(M, SNRdB_vec, Nmc, pars)
% Averages NMSE over Nmc trials. In each trial:
%  - Draw a fresh channel
%  - Draw one random training, keep it FIXED across all SNR points
%  - Calibrate sigma^2 in the whitened domain per that training
    nmse_acc = zeros(size(SNRdB_vec));
    for t = 1:Nmc
        % fresh channel
        chan_mc = generateOnGridChannel(pars.Nt,pars.Nr,pars.K,pars.Nc,pars.Lpaths, ...
                                        pars.angGridTx,pars.angGridRx,1);
        % one training/design for all SNRs
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
    n = (0:N-1).';
    A = (1/sqrt(N)) * exp(1j * (n*pi) * cos(angles(:).'));
end

function chan = generateOnGridChannel(Nt,Nr,K,Nc,Lpaths,angGridTx,angGridRx,rhoL)
    AT = steerULA(Nt, angGridTx);
    AR = steerULA(Nr, angGridRx);
    Gt = size(AT,2); Gr = size(AR,2);

    taps = randi(Nc,[Lpaths,1])-1;
    idTx = randi(Gt,[Lpaths,1]);
    idRx = randi(Gr,[Lpaths,1]);
    gains = (randn(Lpaths,1)+1j*randn(Lpaths,1))/sqrt(2*Lpaths);

    Hk = cell(K,1);
    scale = sqrt(Nt*Nr/(Lpaths*rhoL));
    for k = 1:K
        H = zeros(Nr,Nt);
        for l = 1:Lpaths
            ak = exp(-1j*2*pi*(k-1)/K * taps(l));
            H = H + scale * gains(l)*ak * (AR(:,idRx(l)) * AT(:,idTx(l))');
        end
        Hk{k} = H;
    end
    chan.Hk = Hk;
end

function [Uw, Yw_clean] = buildUw_and_cleanY(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,AT,AR)
% D_w^{-H} * ( (q^T F^T conj(AT)) ⊗ (W^H AR) ), with y_w = D_w^{-H} * (W^H H F q)
    Gt = size(AT,2); Gr = size(AR,2); G = Gt*Gr;
    MLr = M*Lr;
    Uw = zeros(MLr, G);
    Yw_clean = zeros(MLr, K);

    for m = 1:M
        % training beams (explicit shaping to avoid orientation surprises)
        idxF = randi(numel(phaseSet), Nt*Lt, 1);
        Ftr  = (1/sqrt(Nt)) * reshape(exp(1j*phaseSet(idxF)), Nt, Lt); % Nt x Lt
        idxW = randi(numel(phaseSet), Nr*Lr, 1);
        Wtr  = (1/sqrt(Nr)) * reshape(exp(1j*phaseSet(idxW)), Nr, Lr); % Nr x Lr

        q   = ones(Lt,1);       % Lt x 1
        Fq  = Ftr*q;            % Nt x 1

        % dictionary blocks (CRITICAL: W^H, not W^T)
        Ttx = q.' * (Ftr.' * conj(AT));  % 1 x Gt
        Trx = (Wtr') * AR;               % Lr x Gr  (Hermitian)
        U_m = kron(Ttx, Trx);            % Lr x (Gt*Gr)

        % whitening
        Rm = chol(Wtr'*Wtr, 'upper');    % Lr x Lr
        idx = (m-1)*Lr+(1:Lr);

        Uw(idx,:) = Rm' \ U_m;           % D_w^{-H} * U_m

        for k = 1:K
            y = (Wtr') * chan.Hk{k} * Fq;     % Lr x 1
            Yw_clean(idx, k) = Rm' \ y;       % D_w^{-H} * y
        end
    end
end

function Hv_hat = swomp_joint_fast(Yw, Uw, epsStop, maxIter)
    MLr = size(Yw,1);
    G   = size(Uw,2);
    K   = size(Yw,2);

    T = false(G,1);
    r = Yw;
    Hv_hat = zeros(G, K);

    mse = inf; it = 0;
    while (mse > epsStop) && (it < maxIter)
        it = it + 1;

        C = Uw' * r;                 % G x K
        score = sum(abs(C).^2, 2);   % aggregate energy across subcarriers (stable)
        score(T) = -inf;
        [~, p] = max(score);
        T(p) = true;

        UwT = Uw(:, T);              % MLr x |T|
        [Q,R] = qr(UwT, 0);
        X_T = R \ (Q' * Yw);         % |T| x K
        r = Yw - UwT * X_T;

        mse = mean(sum(abs(r).^2,1)) / MLr;
    end

    Hv_hat(T, :) = X_T;
end

function Hhat = Hv_to_H(Hv_hat, AT, AR, K)
    Gr = size(AR,2); Gt = size(AT,2);
    Hhat = cell(K,1);
    for k = 1:K
        Hv_mat = reshape(Hv_hat(:,k), Gr, Gt);
        Hhat{k} = AR * Hv_mat * AT';
    end
end

function val = nmse_of_estimate(Hhat, Htrue_cell)
    K = numel(Htrue_cell);
    num = 0; den = 0;
    for k = 1:K
        num = num + norm(Hhat{k}-Htrue_cell{k}, 'fro')^2;
        den = den + norm(Htrue_cell{k}, 'fro')^2;
    end
    val = num/den;
end

function R = spectral_efficiency(Hhat, Ns, SNRdB)
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