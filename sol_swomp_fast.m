% proj3_swomp_fast_v3.m
% Project 3: SW-OMP – fast, spec-faithful, with SNR calibrated in the whitened domain.
% Plots NMSE vs SNR (M=80,120), NMSE vs M at SNR={-10,-5,0}, and SE vs SNR (M=60).

function sol_swomp_fast()
clc; close all; rng(1);

%% ------------------------- Parameters ------------------------------------
Nt = 32; Nr = 32;
Lt = 1;  Lr = 4;
K  = 16;
Nc = 4;  Lpaths = 4;
Gt = 64; Gr = 64;
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

figure('Name','NMSE vs SNR'); tiledlayout(1,1); ax1=nexttile; hold on; grid on;
for Mi = 1:numel(M_list)
    M = M_list(Mi);
    nmse_curve = nmse_vs_snr_avg(M, SNRdB_vec, Nmc, pars);
    plot(ax1, SNRdB_vec, 10*log10(nmse_curve), 'LineWidth', 1.8, ...
         'DisplayName', sprintf('M=%d',M));
end
xlabel('SNR (dB)'); ylabel('NMSE (dB)');
title('NMSE vs SNR (on-grid AoA/AoD)'); legend('Location','southwest');


%% ------------------------- Fig 2: NMSE vs M ------------------------------
M_sweep = 20:5:100;
SNRdB_set = [-10, -5, 0];

figure('Name','NMSE vs M'); tiledlayout(1,1); ax2=nexttile; hold on; grid on;
for is = 1:numel(SNRdB_set)
    SNRdB = SNRdB_set(is);
    nmse_vsM = zeros(size(M_sweep));

    for im = 1:numel(M_sweep)
        M = M_sweep(im);
        [Uw, Yw_clean] = buildUw_and_cleanY(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,AT,AR);
        sigpow = mean(abs(Yw_clean(:)).^2);
        sigma2 = sigpow / (10^(SNRdB/10));
        Yw = Yw_clean + sqrt(sigma2/2).*(randn(size(Yw_clean))+1j*randn(size(Yw_clean)));

        epsStop = sigma2; 
        maxIter = Lpaths; % change to Lpath
        Hv_hat = swomp_joint_fast(Yw, Uw, epsStop, maxIter);
        Hhat = Hv_to_H(Hv_hat, AT, AR, K);
        nmse_vsM(im) = nmse_of_estimate(Hhat, chan.Hk);
    end

    plot(ax2, M_sweep, 10*log10(nmse_vsM), 'LineWidth',1.8, ...
         'DisplayName',sprintf('SNR=%d dB',SNRdB));
end
xlabel('Training frames, M'); ylabel('NMSE (dB)');
title('NMSE vs number of frames (on-grid AoA/AoD)'); legend('Location','northeast');

%% ------------------------- Fig 3: SE vs SNR (M=60) -----------------------
M = 60;
SNRdB_vec_SE = -15:2.5:10;
SE_curve = zeros(size(SNRdB_vec_SE));
[Uw_SE, Yw_clean_SE] = buildUw_and_cleanY(Nt,Nr,Lt,Lr,M,K,phaseSet,chan,AT,AR);
sigpow_SE = mean(abs(Yw_clean_SE(:)).^2);

for is = 1:numel(SNRdB_vec_SE)
    SNRdB = SNRdB_vec_SE(is);
    sigma2 = sigpow_SE / (10^(SNRdB/10));
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
        % in swomp_joint_fast:
        C = Uw' * r;               % G x K
        score = sum(abs(C).^2, 2); % <-- energy (often more stable)
        score(T) = -inf;
        [~, p] = max(score);
        T(p) = true;

        UwT = Uw(:, T);                  % MLr x |T|
        [Q,R] = qr(UwT, 0);
        X_T = R \ (Q' * Yw);             % |T| x K
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
