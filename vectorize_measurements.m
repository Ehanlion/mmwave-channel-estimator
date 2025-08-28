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

% Stack measurements
y_blocks = cell(M,K);
for m=1:M
    Wm = W_RF(:,:,m)*W_BB(:,:,m);  % [Nr x Ns]
    Fm = F_RF(:,:,m)*F_BB(:,:,m);  % [Nt x Ns]
    for k=1:K
        Ymk = Wm' * H(:,:,k) * Fm * S;   % [Ns x Ns]
        y_blocks{m,k} = Ymk(:);          % vectorize
    end
end
y = cat(1, y_blocks{:});                % [(M*K*Ns*Ns) x 1]

% Add noise to reach target SNR over measurement vector power
sigPow = mean(abs(y).^2);
sigma2 = sigPow / max(10.^(SNRdB/10), eps);
n = sqrt(sigma2/2) * (randn(size(y))+1j*randn(size(y)));
y = y + n;

% Placeholder sensing matrix for SW-OMP (to be constructed by you)
A = [];  % TODO: build dictionary for on-grid sparse recovery (e.g., Kronecker of AoA/AoD ⊗ delay)

% Meta for downstream reconstruction
meta = struct('H_true',H,'W_RF',{W_RF},'W_BB',{W_BB},'F_RF',{F_RF},'F_BB',{F_BB},...
              'S',S,'sigma2',sigma2,'dims',[Nr Nt K Ns M]);
end
