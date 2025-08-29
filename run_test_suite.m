% run_test_suite.m - A script to test and debug the mmWave project functions.
%
% Description:
%   This script runs a single instance of the channel estimation pipeline to
%   verify that all helper functions are working correctly. It prints the
%   dimensions of key matrices and the final NMSE to the command window,
%   allowing for easy debugging and validation of each step.

% --- SETUP ---
clear; clc; close all;
% Add the functions directory to the MATLAB path
addpath('functions');

fprintf('--- Starting Project Function Test Suite ---\n\n');

% --- 1. Test initialize_parameters ---
fprintf('1. Testing initialize_parameters()...\n');
params = initialize_parameters();
if isstruct(params) && isfield(params, 'Nt')
    fprintf('   SUCCESS: Parameters loaded.\n\n');
else
    fprintf('   FAILURE: Could not load parameters.\n');
    return;
end

% --- 2. Test generate_mmwave_channel ---
fprintf('2. Testing generate_mmwave_channel()...\n');
[H_true_freq, ~, ~, ~, ~, ~, ~] = generate_mmwave_channel(params);
if iscell(H_true_freq) && all(size(H_true_freq{1}) == [params.Nr, params.Nt])
    fprintf('   SUCCESS: Channel generated.\n');
    fprintf('   - H_true_freq size: %d x 1 cell array\n', length(H_true_freq));
    fprintf('   - Matrix size inside cell: %d x %d\n\n', size(H_true_freq{1}, 1), size(H_true_freq{1}, 2));
else
    fprintf('   FAILURE: Channel generation failed.\n');
    return;
end

% --- 3. Setup for Vectorization and SW-OMP ---
fprintf('3. Setting up training signals and noise...\n');
M = 80; % Number of training frames for this test
SNR_dB = 0; % Test SNR in dB
snr_linear = 10^(SNR_dB / 10);
noise_var = 1 / snr_linear;

% Generate random training precoders and combiners with quantized phases
F_train = cell(M, 1);
W_train = cell(M, 1);
q_train = cell(M, 1);
phase_set = (0:(2^params.Nq - 1)) * 2 * pi / (2^params.Nq);

for m = 1:M
    % Precoders
    rand_phases_f = phase_set(randi(length(phase_set), params.Nt, params.Lt));
    F_train{m} = (1/sqrt(params.Nt)) * exp(1j * rand_phases_f);
    % Combiners
    rand_phases_w = phase_set(randi(length(phase_set), params.Nr, params.Lr));
    W_train{m} = (1/sqrt(params.Nr)) * exp(1j * rand_phases_w);
    % Symbols
    q_train{m} = (randn(params.Lt, 1) + 1j*randn(params.Lt, 1))/sqrt(2);
end
t_train = ones(params.K, M); % Simple pilot symbols for testing
fprintf('   SUCCESS: Training signals created.\n\n');

% --- 4. Test vectorize_measurements ---
fprintf('4. Testing vectorize_measurements()...\n');
[y_vec, Phi, Psi, h_v_true] = vectorize_measurements(params, H_true_freq, F_train, W_train, q_train, t_train, noise_var);
fprintf('   SUCCESS: Vectorization complete.\n');
fprintf('   - Phi matrix size: %d x %d\n', size(Phi, 1), size(Phi, 2));
fprintf('   - Psi matrix size: %d x %d\n', size(Psi, 1), size(Psi, 2));
fprintf('   - y_vec{1} vector size: %d x %d\n\n', size(y_vec{1}, 1), size(y_vec{1}, 2));

% --- 5. Test swomp ---
fprintf('5. Testing swomp()...\n');
[h_v_est, T_est] = swomp(params, y_vec, Phi, Psi, noise_var);
fprintf('   SUCCESS: SW-OMP algorithm completed.\n');
fprintf('   - Estimated support size: %d paths\n', length(T_est));
fprintf('   - Estimated support indices: %s\n\n', mat2str(T_est));

% --- 6. Reconstruct Channel and Test calculate_nmse ---
fprintf('6. Reconstructing channel and testing calculate_nmse()...\n');
H_est_freq = cell(params.K, 1);
for k = 1:params.K
    % Reconstruct the channel matrix from the sparse vector
    vec_H_est = Psi * h_v_est{k};
    H_est_freq{k} = reshape(vec_H_est, params.Nr, params.Nt);
end

nmse_db = calculate_nmse(H_true_freq, H_est_freq);
fprintf('   SUCCESS: NMSE calculated.\n');
fprintf('   - Final NMSE: %.2f dB\n\n', nmse_db);

fprintf('--- Test Suite Finished Successfully ---\n');
