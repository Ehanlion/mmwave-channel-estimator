% channel_params.m
% Create and return a params struct with project-wide defaults.
% Usage:
%   p = channel_params();                     % defaults
%   p = channel_params("M",120,"SNRdB",0);    % overrides (name-value)
function params = channel_params(varargin)
params = struct();

% --- Arrays & streams ---
params.Nt    = 32; params.Nr = 16;
params.Nrf_t = 4;  params.Nrf_r = 4;
params.Ns    = 2;

% --- OFDM / bandwidth ---
params.K     = 64;
params.fc    = 28e9;
params.BW    = 400e6;
params.fs    = params.BW;

% --- Channel sparsity (on-grid) ---
params.L        = 3;
params.Ng_tx    = 64;
params.Ng_rx    = 64;
params.Ng_tau   = 16;
params.angles_on_grid = true;

% --- Training ---
params.M        = 80;
params.pilotPow = 1;

% --- Noise / SNR ---
params.SNRdB_list = -15:5:10;
params.SNRdB      = 0;

% --- MC ---
params.Nmc = 50;

% --- Debug / logging ---
params.verbose     = false;
params.debug_dump  = false;
params.debug_dir   = fullfile(pwd,'debug_dumps');
params.log_dir     = fullfile(pwd,'logs');
params.random_seed = 233;

% --- Reconstruction ---
params.mock_reconstruction = false;   % <--- use real SW-OMP now
params.max_iters = params.L;          % SW-OMP iterations cap
params.stop_tol  = 1e-3;              % relative residual improvement threshold

% overrides
if ~isempty(varargin)
    assert(mod(numel(varargin),2)==0,"Use name-value pairs.");
    for i=1:2:numel(varargin)
        params.(string(varargin{i})) = varargin{i+1};
    end
end
end