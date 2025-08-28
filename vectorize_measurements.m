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
Nr = params.Nr; Nt = params.Nt; K = params.K; M = params.M; Ns = params.Ns;
pilotPow = params.pilotPow;

% Pilot S (Ns x Ns): identity per shot (one shot)
S = sqrt(pilotPow)*eye(Ns);

% Stack measurements (signal only)
y_blocks = cell(M,K);
for m=1:M
    Wm = W_RF(:,:,m)*W_BB(:,:,m);  % [Nr x Ns]
    Fm = F_RF(:,:,m)*F_BB(:,:,m);  % [Nt x Ns]
    for k=1:K
        Ymk = Wm' * H(:,:,k) * Fm * S;   % [Ns x Ns]
        y_blocks{m,k} = Ymk(:);
    end
end
y_clean = cat(1, y_blocks{:});            % signal-only measurements
sigPow = mean(abs(y_clean).^2);

% Add noise to reach target SNR over measurement vector power
sigma2 = sigPow / max(10.^(SNRdB/10), eps);
n = sqrt(sigma2/2) * (randn(size(y_clean))+1j*randn(size(y_clean)));
y = y_clean + n;

% Achieved SNR sanity (due to finite sample randomness)
achSNR = 10*log10(mean(abs(y_clean).^2)/mean(abs(n).^2));

% Placeholder sensing matrix for SW-OMP (to be constructed)
A = [];

% Debug prints
dbg(params,'[vectorize] dims: Nr=%d Nt=%d Ns=%d K=%d M=%d', Nr,Nt,Ns,K,M);
dbg(params,'[vectorize] |y|^2 mean (signal)=%.3e, sigma^2 set=%.3e, target SNR=%.2f dB, achieved SNR≈%.2f dB', ...
    sigPow, sigma2, SNRdB, achSNR);
dbg(params,'[vectorize] y length: %d (M*K*Ns^2=%d)', numel(y), M*K*Ns*Ns);

% Meta for downstream reconstruction & optional dump
meta = struct('H_true',H,'W_RF',{W_RF},'W_BB',{W_BB},'F_RF',{F_RF},'F_BB',{F_BB},...
              'S',S,'sigma2',sigma2,'dims',[Nr Nt K Ns M],'achSNR',achSNR);

if isfield(params,'debug_dump') && params.debug_dump
    if ~exist(params.debug_dir,'dir'); mkdir(params.debug_dir); end
    stamp = char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
    save(fullfile(params.debug_dir, sprintf('vecdump_%s.mat',stamp)), ...
         'H','W_RF','W_BB','F_RF','F_BB','y','y_clean','n','meta','params','SNRdB','-v7.3');
    dbg(params,'[vectorize] dump saved.');
end
end