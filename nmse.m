% nmse.m
% NMSE across subcarriers: sum_k ||H - Hhat||_F^2 / sum_k ||H||_F^2
function [val, perK] = nmse(H_true, H_hat, params)
if nargin < 3, params = struct(); end
assert(ndims(H_true)==3 && all(size(H_true)==size(H_hat)), 'Size mismatch.');

K = size(H_true,3);
num = 0; den = 0;
perK = zeros(1,K);
for k=1:K
    ek = norm(H_true(:,:,k)-H_hat(:,:,k),'fro')^2;
    hk = norm(H_true(:,:,k),'fro')^2;
    perK(k) = ek / max(hk,eps);
    num = num + ek;
    den = den + hk;
end
val = num / max(den, eps);

if isfield(params,'verbose') && params.verbose
    fprintf('[nmse] NMSE=%.3e | mean(perK)=%.3e, std(perK)=%.3e\n', val, mean(perK), std(perK));
end
end