% SOLUTION FILE, plotting Spectral Efficiency vs. SNR
% spectral_vs_snr.m
% Experiment: Spectral Efficiency vs SNR at fixed M=60. Produces a plot.
clear; close all; clc;

params = channel_params("M",60, "mock_reconstruction", true); % set false after SW-OMP done
params.Nmc = 20;
SNRdB_list = -15:5:10;

rates = zeros(size(SNRdB_list));

for iS = 1:numel(SNRdB_list)
    snrdb = SNRdB_list(iS);
    acc = 0;
    for t=1:params.Nmc
        [H,~] = gen_channel(params);
        [W_RF,W_BB,F_RF,F_BB] = build_training(params);
        [y, A, meta] = vectorize_measurements(H, W_RF,W_BB,F_RF,F_BB, params, snrdb); %#ok<ASGLU>
        [rec,~] = swomp(y, A, meta, params);
        acc = acc + spectral_efficiency(rec.H_hat, W_RF,W_BB,F_RF,F_BB, params, snrdb);
    end
    rates(iS) = acc/params.Nmc;
end

figure; grid on; hold on;
plot(SNRdB_list, rates, '-o');
xlabel('SNR [dB]'); ylabel('Spectral Efficiency [bit/s/Hz]');
title('Spectral Efficiency vs SNR (M=60)');
