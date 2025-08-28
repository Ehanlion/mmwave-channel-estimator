% SOLUTION FILE for NMSE vs. SNR
% exp_nmse_vs_snr.m
% Experiment: NMSE vs SNR for M in {80, 120}. Produces a plot.
clear; close all; clc;

% --- Logging setup ---
params0 = channel_params();
if ~exist(params0.log_dir,'dir'); mkdir(params0.log_dir); end
logfile = fullfile(params0.log_dir, sprintf('exp_nmse_vs_snr_%s.log', char(datetime('now','Format','yyyyMMdd_HHmmss'))));
diary(logfile); diary on;

fprintf('== EXP: NMSE vs SNR ==\n');

% --- Base params & sweeps ---
M_list = [80 120];
SNRdB_list = -15:5:10;

nmse_curves = zeros(numel(M_list), numel(SNRdB_list));

for iM = 1:numel(M_list)
    params = channel_params("M", M_list(iM), ...
                            "random_seed", 100+iM, ...
                            "mock_reconstruction", true, ...   % set false after SW-OMP done
                            "verbose", true, ...
                            "debug_dump", false);
    fprintf('--- M=%d ---\n', params.M);
    for iS = 1:numel(SNRdB_list)
        snrdb = SNRdB_list(iS);
        acc = 0;
        for t=1:params.Nmc
            [H,~] = gen_channel(params);
            [W_RF,W_BB,F_RF,F_BB] = build_training(params);
            [y, A, meta] = vectorize_measurements(H, W_RF,W_BB,F_RF,F_BB, params, snrdb); %#ok<ASGLU>
            [rec,info] = swomp(y, A, meta, params); %#ok<NASGU>
            acc = acc + nmse(H, rec.H_hat, params);
        end
        nmse_curves(iM,iS) = acc/params.Nmc;
        fprintf('M=%d | SNR=%+3.0f dB -> NMSE=%.3e\n', params.M, snrdb, nmse_curves(iM,iS));
    end
end

% Flat-line warning (likely MOCK mode)
rngvar = std(nmse_curves(:));
if max(abs(nmse_curves(:))) < 1e-6 || rngvar < 1e-8
    warning('NMSE curve is ~flat (≈0). MOCK reconstruction is probably enabled (params.mock_reconstruction=true).');
end

% Plot
figure; hold on; grid on;
plot(SNRdB_list, nmse_curves(1,:), '-o', 'DisplayName', 'M=80');
plot(SNRdB_list, nmse_curves(2,:), '-s', 'DisplayName', 'M=120');
xlabel('SNR [dB]'); ylabel('NMSE');
title('NMSE vs SNR (on-grid angles)');
legend('Location','southwest');

% Save results
outdir = fullfile(pwd,'fig'); if ~exist(outdir,'dir'); mkdir(outdir); end
saveas(gcf, fullfile(outdir,'nmse_vs_snr_M80_M120.png'));
save(fullfile(outdir,'nmse_vs_snr_data.mat'), 'SNRdB_list','M_list','nmse_curves');

fprintf('Logs written to: %s\n', logfile);
diary off;