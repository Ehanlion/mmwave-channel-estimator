% SOLUTION FILE for NMSE vs. SNR
% exp_nmse_vs_snr.m
% Experiment: NMSE vs SNR for M in {80, 120}. Produces a plot.
clear; close all; clc;

% Base params
base = channel_params();
base.Nmc = 20;                 % keep light; increase later
M_list = [80 120];
SNRdB_list = -15:5:10;

nmse_curves = zeros(numel(M_list), numel(SNRdB_list));

for iM = 1:numel(M_list)
    params = channel_params("M", M_list(iM), "random_seed", 100+iM, ...
                            "mock_reconstruction", true);  % set false after SW-OMP done
    for iS = 1:numel(SNRdB_list)
        snrdb = SNRdB_list(iS);
        acc = 0;
        for t=1:params.Nmc
            [H,~] = gen_channel(params);
            [W_RF,W_BB,F_RF,F_BB] = build_training(params);
            [y, A, meta] = vectorize_measurements(H, W_RF,W_BB,F_RF,F_BB, params, snrdb); %#ok<ASGLU>
            [rec,~] = swomp(y, A, meta, params);
            acc = acc + nmse(H, rec.H_hat);
        end
        nmse_curves(iM,iS) = acc/params.Nmc;
    end
end

figure; hold on; grid on;
plot(SNRdB_list, nmse_curves(1,:), '-o', 'DisplayName', 'M=80');
plot(SNRdB_list, nmse_curves(2,:), '-s', 'DisplayName', 'M=120');
xlabel('SNR [dB]'); ylabel('NMSE');
title('NMSE vs SNR (on-grid angles)');
legend('Location','southwest');
