% vectorize_measurements.m
% Form stacked training measurements y across K subcarriers and M frames.
% This skeleton returns:
%   y   : stacked vector of all measurements (complex)
%   A   : sensing matrix (TO FILL for SW-OMP) -> placeholder empty in skeleton
%   meta: struct carrying H_true and training for downstream reconstruction
%
% Measurement model per frame m, subcarrier k:
%   Y_mk = W_m^H * H[k] * F_m * S,  with S (Ns x Ns_pilot). Here we use S = I_Ns, Ns_pilot = Ns,
%   then we stack vec(Y_mk) over (k,m).

function [y, A, meta] = vectorize_measurements(H, W_RF,W_BB,F_RF,F_BB, params, SNRdB)
Nr=params.Nr; Nt=params.Nt; K=params.K; M=params.M; Ns=params.Ns;

% Pilot (kept identity for clarity)
S = eye(Ns);

y_blocks = cell(M,K);          % signal-only (after combining)
n_blocks = cell(M,K);          % noise after combining (colored by W)

for m=1:M
    Wm = W_RF(:,:,m)*W_BB(:,:,m);    % [Nr x Ns]
    Fm = F_RF(:,:,m)*F_BB(:,:,m);    % [Nt x Ns]
    for k=1:K
        Ymk_sig = Wm' * H(:,:,k) * Fm * S;             % [Ns x Ns]
        y_blocks{m,k} = Ymk_sig(:);

        % Receiver noise BEFORE combining, then combine:
        Nmk = (randn(Nr,Ns)+1j*randn(Nr,Ns))/sqrt(2);  % unit-variance
        Ymk_n = Wm' * Nmk;                              % [Ns x Ns]
        n_blocks{m,k} = Ymk_n(:);
    end
end

y_sig = cat(1,y_blocks{:});
n_raw = cat(1,n_blocks{:});

% Scale noise to hit target SNR over stacked measurements
Ps  = mean(abs(y_sig).^2);
Pn0 = mean(abs(n_raw).^2);
snr_lin = 10^(SNRdB/10);
alpha = sqrt( (Ps/(snr_lin*Pn0)) );      % scale n_raw -> target noise power
n = alpha * n_raw;
y = y_sig + n;

% Diagnostics
achSNR = 10*log10(mean(abs(y_sig).^2)/mean(abs(n).^2));

A = []; % dictionary built implicitly in swomp

dbg(params,'[vectorize] dims: Nr=%d Nt=%d Ns=%d K=%d M=%d', Nr,Nt,Ns,K,M);
dbg(params,'[vectorize] |y_sig|^2 mean=%.3e, |n|^2 mean=%.3e, target=%.2f dB, achieved≈%.2f dB', ...
    Ps, mean(abs(n).^2), SNRdB, achSNR);
dbg(params,'[vectorize] y length: %d (M*K*Ns^2=%d)', numel(y), M*K*Ns*Ns);

meta = struct('H_true',H,'W_RF',{W_RF},'W_BB',{W_BB},'F_RF',{F_RF},'F_BB',{F_BB},...
              'S',S,'dims',[Nr Nt K Ns M],'achSNR',achSNR);
end
