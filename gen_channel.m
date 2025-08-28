% gen_channel.m
% Generate a wideband mmWave MIMO channel H(:,:,k) with L on-grid paths.
% Returns:
%   H : [Nr x Nt x K] frequency-domain channel over K subcarriers
%   chinfo : struct with AoAs/AoDs/delays/gains
function [H, chinfo] = gen_channel(params)
seed_rng(params.random_seed); %#ok<*NASGU>

Nt = params.Nt; Nr = params.Nr; K = params.K; L = params.L;
Ng_tx = params.Ng_tx; Ng_rx = params.Ng_rx; Ng_tau = params.Ng_tau;
fs  = params.fs;

% Uniform linear arrays (half-wavelength spacing)
ula_steer = @(N,theta) exp(1j*pi*(0:N-1).'*sin(theta)) / sqrt(N);

% On-grid angles & delays
tx_grid = asin(linspace(-1,1,Ng_tx));  % AoD grid (radians)
rx_grid = asin(linspace(-1,1,Ng_rx));  % AoA grid (radians)
tau_grid = (0:Ng_tau-1).'/(fs);        % delay grid (seconds), simple uniform grid

% Random L-path support on grids
tx_idx = randi(Ng_tx,L,1);
rx_idx = randi(Ng_rx,L,1);
tau_idx= randi(Ng_tau,L,1);

aod = tx_grid(tx_idx);
aoa = rx_grid(rx_idx);
tau = tau_grid(tau_idx);

% Complex path gains (normalized)
alpha = (randn(L,1)+1j*randn(L,1))/sqrt(2*L);

% Subcarrier frequencies
f = linspace(-params.BW/2, params.BW/2, K);  % baseband subcarrier freqs

H = zeros(Nr,Nt,K);
for k = 1:K
    Hk = zeros(Nr,Nt);
    for ell = 1:L
        ar = ula_steer(Nr, aoa(ell));
        at = ula_steer(Nt, aod(ell));
        freq_phase = exp(-1j*2*pi*f(k)*tau(ell));
        Hk = Hk + alpha(ell) * freq_phase * (ar * at');
    end
    H(:,:,k) = Hk;
end

% Normalize average Frobenius norm across K to 1 (optional)
pow = mean(arrayfun(@(kk) norm(H(:,:,kk),'fro')^2, 1:K));
if pow > 0
    H = H / sqrt(pow);
end

chinfo = struct('aoa',aoa,'aod',aod,'tau',tau,'alpha',alpha,...
                'tx_idx',tx_idx,'rx_idx',rx_idx,'tau_idx',tau_idx);
end
