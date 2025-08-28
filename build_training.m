% build_training.m
% Create random hybrid training precoders/combiners for M frames.
% Returns 4 arrays:
%   W_RF: [Nr x Nrf_r x M],  W_BB: [Nrf_r x Ns x M]
%   F_RF: [Nt x Nrf_t x M],  F_BB: [Nrf_t x Ns x M]
function [W_RF,W_BB,F_RF,F_BB] = build_training(params)
seed_rng(params.random_seed);

Nt=params.Nt; Nr=params.Nr; M=params.M;
Nrf_t=params.Nrf_t; Nrf_r=params.Nrf_r; Ns=params.Ns;

W_RF = zeros(Nr,Nrf_r,M);   F_RF = zeros(Nt,Nrf_t,M);
W_BB = zeros(Nrf_r,Ns,M);   F_BB = zeros(Nrf_t,Ns,M);

for m=1:M
    W_RF(:,:,m) = exp(1j*2*pi*rand(Nr,Nrf_r))/sqrt(Nr);
    F_RF(:,:,m) = exp(1j*2*pi*rand(Nt,Nrf_t))/sqrt(Nt);

    tmp = (randn(Nrf_t,Ns)+1j*randn(Nrf_t,Ns))/sqrt(2*Nrf_t);
    F_BB(:,:,m) = tmp;
    scale = sqrt(trace((F_RF(:,:,m)*F_BB(:,:,m))'*(F_RF(:,:,m)*F_BB(:,:,m)))/Ns);
    F_BB(:,:,m) = F_BB(:,:,m)/max(scale,eps);

    [Q,~] = qr((randn(Nrf_r,Nrf_r)+1j*randn(Nrf_r,Nrf_r))/sqrt(2),0);
    W_BB(:,:,m) = Q(:,1:Ns);
end

dbg(params,'[build_training] Nr=%d Nt=%d Nrf_r=%d Nrf_t=%d Ns=%d M=%d', ...
    Nr,Nt,Nrf_r,Nrf_t,Ns,M);
end
