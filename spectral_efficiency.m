% spectral_efficiency.m
% Compute average per-subcarrier achievable spectral efficiency (bit/s/Hz)
% using hybrid precoder/combiner derived from H_hat (skeleton uses training beams).
% Inputs:
%   H_hat, W_RF/W_BB, F_RF/F_BB, params, snr_db
% Outputs:
%   R (scalar): average over K
function R = spectral_efficiency(H_hat, W_RF,W_BB,F_RF,F_BB, params, snr_db)
K  = size(H_hat,3); Ns = params.Ns;
rho = 10.^(snr_db/10);   % SNR (linear) per Rx noise power = 1

% Skeleton: use training beams as effective combiner/precoder during rate eval.
% Replace with your data-mode hybrid design (e.g., SVD-based BB with RF constraints).
sumRate = 0;
for k=1:K
    % Effective baseband channel Heff = W^H * H * F
    % Aggregate an average across M training frames (skeleton simplification).
    Heff_accum = zeros(Ns,Ns);
    for m=1:params.M
        Wm = W_RF(:,:,m)*W_BB(:,:,m);
        Fm = F_RF(:,:,m)*F_BB(:,:,m);
        Heff_accum = Heff_accum + (Wm')*H_hat(:,:,k)*Fm;
    end
    Heff = Heff_accum / params.M;

    % Achievable rate with equal power and no interference cancelation:
    % Rk = log2 det( I + (rho/Ns) * Heff*Heff^H )
    G = Heff*Heff';
    Rk = real(log2(det( eye(Ns) + (rho/max(Ns,1))*G )));
    sumRate = sumRate + max(Rk,0);
end
R = sumRate / K;
end
