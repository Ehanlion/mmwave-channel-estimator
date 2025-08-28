% nmse.m
% NMSE across subcarriers: sum_k ||H - Hhat||_F^2 / sum_k ||H||_F^2
function val = nmse(H_true, H_hat)
assert(ndims(H_true)==3 && all(size(H_true)==size(H_hat)), 'Size mismatch.');

K = size(H_true,3);
num = 0; den = 0;
for k=1:K
    num = num + norm(H_true(:,:,k)-H_hat(:,:,k),'fro')^2;
    den = den + norm(H_true(:,:,k),'fro')^2;
end
val = num / max(den, eps);
end
