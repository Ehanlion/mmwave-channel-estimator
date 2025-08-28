% channel_params.m
% Create and return a params struct with project-wide defaults.
% Usage:
%   p = channel_params();                     % defaults
%   p = channel_params("M",120,"SNRdB",0);    % overrides (name-value)
function params = channel_params(varargin)
params = struct();

% --- Arrays & streams ---
params.Nt    = 32;     % TX antennas
params.Nr    = 16;     % RX antennas
params.Nrf_t = 4;      % TX RF chains
params.Nrf_r = 4;      % RX RF chains
params.Ns    = 2;      % data streams (<= min(Nrf_t, Nrf_r))

% --- OFDM / bandwidth ---
params.K     = 64;     % subcarriers
params.fc    = 28e9;   % carrier [Hz]
params.BW    = 400e6;  % bandwidth [Hz]
params.fs    = params.BW;  % sample rate (simple wideband model)

% --- Channel sparsity (on-grid) ---
params.L        = 3;               % paths
params.Ng_tx    = 64;              % AoD grid size
params.Ng_rx    = 64;              % AoA grid size
params.Ng_tau   = 16;              % delay grid size
params.angles_on_grid = true;      % on-grid assumption (required)

% --- Training (to be swept in experiments) ---
params.M        = 80;              % training frames (also used as SW-OMP iters bound)
params.pilotPow = 1;               % pilot symbol power

% --- Noise / SNR (over subcarrier average) ---
params.SNRdB_list = -15:5:10;
params.SNRdB      = 0;             % default if single run

% --- MC ---
params.Nmc = 50;                   % Monte Carlo trials (tune in experiments)

% --- Flags / skeleton helpers ---
params.mock_reconstruction = true; % if true, SW-OMP returns oracle H_hat = H_true (skeleton)
params.random_seed         = 233;  % reproducibility

% --- Apply name-value overrides ---
if ~isempty(varargin)
    assert(mod(numel(varargin),2)==0,"Use name-value pairs.");
    for i=1:2:numel(varargin)
        params.(string(varargin{i})) = varargin{i+1};
    end
end
end
