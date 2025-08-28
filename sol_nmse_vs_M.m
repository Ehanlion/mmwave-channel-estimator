% SOLUTION FILE, NMSE vs. M
% exp_nmse_vs_M.m
% Experiment: NMSE vs M for SNR in {-10, -5, 0} dB. Produces a plot.
clear; close all; clc;

base = channel_params();
base.Nmc = 20;

M_vals = 20:5:100;
SNRdB_choices = [-10 -5 0];
nmse_mat = zeros(numel(SNRdB_choices), numel(M_vals));

for iS = 1:numel(SNRdB_choices)
    snrdb = SNRdB_choices(iS);
    for iM = 1:numel(M_vals)
        M = M_vals(iM);
        params = channel_params("M", M, "random_seed", 200+iM, ...
                                "mock_reconstruction", true); % set false after SW-OMP done
        acc = 0;
        for t=1:params.Nmc
            [H,~] = gen_channel(params);
            [W_RF,W_BB,F_RF,F_BB] = build_training(params);
            [y, A, meta] = vectorize_measurements(H, W_RF,W_BB,F_RF,F_BB, params, snrdb); 
            [rec,~] = swomp(y, A, meta, params);
            acc = acc + nmse(H, rec.H_hat);
        end
        nmse_mat(iS,iM) = acc/params.Nmc;
    end
end

figure; hold on; grid on;
for iS=1:numel(SNRdB_choices)
    plot(M_vals, nmse_mat(iS,:), '-o', 'DisplayName', sprintf('SNR=%+ddB',SNRdB_choices(iS)));
end
xlabel('Training frames M'); ylabel('NMSE');
title('NMSE vs M at fixed SNRs (on-grid angles)');
legend('Location','northeast');
