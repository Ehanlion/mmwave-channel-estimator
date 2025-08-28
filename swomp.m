% swomp.m
% Simultaneous Weighted OMP (skeleton).
% Inputs:
%   y, A, meta, params
% Outputs:
%   rec : struct with fields
%       .H_hat  [Nr x Nt x K]  (reconstructed channel)
%       .x_hat  (sparse coeffs if using dictionary; optional)
%       .support (selected atom indices; optional)
%   info: struct with diagnostics

function [rec, info] = swomp(y, ~, meta, params)
% SW-OMP with whitening (per Alg.1).
% y     : [K*M*Ns^2 x 1] stacked by k (outer), m (inner)
% meta  : W_RF,W_BB,F_RF,F_BB,dims=[Nr Nt K Ns M]
% params: Ng_rx,Ng_tx,Ng_tau, fs, BW; optional max_iters, stop_tol, verbose

% ---- dims & training ----
Nr=meta.dims(1); Nt=meta.dims(2); K=meta.dims(3); Ns=meta.dims(4); M=meta.dims(5);
W_RF=meta.W_RF; W_BB=meta.W_BB; F_RF=meta.F_RF; F_BB=meta.F_BB;

% ---- grids & steering ----
Ng_rx=params.Ng_rx; Ng_tx=params.Ng_tx; Ng_tau=params.Ng_tau;
fs=params.fs; BW=params.BW;

ula=@(N,th) exp(1j*pi*(0:N-1).'*sin(th))/sqrt(N);
txg=asin(linspace(-1,1,Ng_tx)); rxg=asin(linspace(-1,1,Ng_rx));
taug=(0:Ng_tau-1).'/(fs);
f = linspace(-BW/2,BW/2,K);                % [1xK]
D = exp(-1j*2*pi*(taug*f));                % [Ng_tau x K]

A_tx=zeros(Nt,Ng_tx); for i=1:Ng_tx, A_tx(:,i)=ula(Nt,txg(i)); end
A_rx=zeros(Nr,Ng_rx); for i=1:Ng_rx, A_rx(:,i)=ula(Nr,rxg(i)); end
A_ang = kron(conj(A_tx),A_rx);             % [Nr*Nt x Qang]
Qang  = Ng_rx*Ng_tx;

% ---- build T_stack = blkstack_m(Fm^T ⊗ Wm^H) and whitening Wblk ----
blk = Ns*Ns; Lm = M*blk;
T_stack = zeros(Lm,Nr*Nt);
Wblk    = eye(Lm);
for m=1:M
    Wm = W_RF(:,:,m)*W_BB(:,:,m);        % [Nr x Ns]
    Fm = F_RF(:,:,m)*F_BB(:,:,m);        % [Nt x Ns]
    Tm = kron(Fm.', Wm');                % [Ns^2 x Nr*Nt]
    rows = (m-1)*blk + (1:blk);
    T_stack(rows,:)=Tm;

    % Whitening for colored noise: Cov{vec(Wm' n)} = I ⊗ (Wm'Wm)
    G  = Wm'*Wm;  G = (G+G')/2;          % hermitian safeguard
    Gi = sqrtm(G)\eye(Ns);               % (Wm'Wm)^(-1/2)
    Wm_wh = kron(eye(Ns), Gi);           % [Ns^2 x Ns^2]
    Wblk(rows,rows) = Wm_wh;
end

% Whitened operator
B = (Wblk*T_stack)*A_ang;                % [Lm x Qang]
bnrm = sum(abs(B).^2,1).'; bnrm(bnrm==0)=eps;
bnrm_full = repmat(bnrm,Ng_tau,1);       % [Qang*Ng_tau x 1]

% ---- split & whiten y ----
yK=cell(K,1); rK=cell(K,1);
for k=1:K
    idx=(k-1)*Lm+(1:Lm);
    yK{k}=Wblk*y(idx);                   % whitened measurements
    rK{k}=yK{k};
end
y_stack = cell2mat(yK);

% ---- params ----
maxIters = min(getdef(params,'max_iters',params.L), Ng_tau*Ng_rx*Ng_tx);
stop_tol = getdef(params,'stop_tol',1e-3);
verbose  = isfield(params,'verbose') && params.verbose;

S_lin=[]; S_q=[]; S_l=[]; xS=[];
res_hist=zeros(maxIters+1,1); res_hist(1)=norm(y_stack)^2;

if verbose, fprintf('[SW-OMP] K=%d, M=%d, Ns=%d, Qang=%d, Ng_tau=%d, Lm=%d\n',K,M,Ns,Qang,Ng_tau,Lm); end

% =========================
%      Greedy iterations
% =========================
for t=1:maxIters
    % Distributed correlations per k in whitened domain
    S = zeros(Qang,K);
    for k=1:K
        S(:,k) = B' * rK{k};            % [Qang x 1]
    end
    % Coherent aggregation across subcarriers (delay focusing)
    C = S * conj(D).';                  % [Qang x Ng_tau]
    scr = abs(C(:)) ./ sqrt(bnrm_full); % normalize by dictionary column energy

    if ~isempty(S_lin), scr(S_lin)=-inf; end
    [~,idx_lin]=max(scr);
    l = ceil(idx_lin/Qang);  q=idx_lin-(l-1)*Qang;

    S_lin(end+1,1)=idx_lin; S_q(end+1,1)=q; S_l(end+1,1)=l;
    if verbose, fprintf('[SW-OMP] it=%d pick: q=%d/%d, l=%d/%d, score=%.3e\n',t,q,Qang,l,Ng_tau,scr(idx_lin)); end

    % Build whitened design AS and solve LS jointly
    AS=zeros(K*Lm,t);
    for j=1:t
        aj=A_ang(:,S_q(j));
        lj=S_l(j);
        col=zeros(K*Lm,1);
        for k=1:K
            uk = (Wblk*T_stack) * ( exp(-1j*2*pi*f(k)*taug(lj)) * aj ); % whitened column
            idx=(k-1)*Lm+(1:Lm);
            col(idx)=uk;
        end
        AS(:,j)=col;
    end
    xS = AS \ y_stack;

    % Update residuals in whitened domain
    r_stack = y_stack - AS*xS;
    for k=1:K
        idx=(k-1)*Lm+(1:Lm);
        rK{k}=r_stack(idx);
    end
    res_hist(t+1)=norm(r_stack)^2;

    if res_hist(t)>0
        imp=(res_hist(t)-res_hist(t+1))/res_hist(t);
        if verbose, fprintf('[SW-OMP] it=%d residual=%.3e, Δ=%.3e\n',t,res_hist(t+1),imp); end
        if imp<stop_tol, break; end
    end
end
T=numel(S_lin);

% ---- reconstruct H_hat over K (unwhitened model) ----
H_hat=zeros(Nr,Nt,K);
for k=1:K
    hvec=zeros(Nr*Nt,1);
    for j=1:T
        hvec = hvec + (exp(-1j*2*pi*f(k)*taug(S_l(j))) * A_ang(:,S_q(j))) * xS(j);
    end
    H_hat(:,:,k)=reshape(hvec,Nr,Nt);
end

rec.H_hat=H_hat; rec.x_hat=xS; rec.support=[S_q,S_l];
info.residual_hist=res_hist(1:min(T+1,end));
info.Qang=Qang; info.maxIters=maxIters; info.stop_tol=stop_tol;
end

function v = getdef(s, field, dflt)
    if isfield(s,field) && ~isempty(s.(field)), v = s.(field); else, v = dflt; end
end
