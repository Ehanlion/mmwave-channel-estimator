function [y_vec, Phi, Psi, h_v] = vectorize_measurements(params, H_freq, F_train, W_train, q_train, t_train, noise_var)
    % vectorize_measurements - Formulates the compressive sensing problem y = Phi * Psi * h_v.
    %
    % Syntax: [y_vec, Phi, Psi, h_v] = vectorize_measurements(params, H_freq, F_train, W_train, q_train, t_train, noise_var)
    %
    % Inputs:
    %   params    - Struct with simulation parameters.
    %   H_freq    - K x 1 cell array of frequency-domain channel matrices.
    %   F_train   - M x 1 cell array of training precoders.
    %   W_train   - M x 1 cell array of training combiners.
    %   q_train   - M x 1 cell array of frequency-flat symbol vectors.
    %   t_train   - K x M matrix of pilot symbols.
    %   noise_var - Variance of the AWGN.
    %
    % Outputs:
    %   y_vec - K x 1 cell array of stacked measurement vectors for each subcarrier.
    %   Phi   - (M*Lr) x (Nt*Nr) sensing matrix.
    %   Psi   - (Nt*Nr) x (Gt*Gr) dictionary matrix.
    %   h_v   - K x 1 cell array of true sparse channel vectors.
    %
    % Description:
    %   This function implements the vectorization process from [R-4, Sec. III-A].
    %   It constructs the overall sensing matrix (Phi), the dictionary matrix (Psi),
    %   the received measurement vectors (y_vec), and the true sparse channel
    %   vectors (h_v) needed for the SW-OMP algorithm.

    % Extract parameters
    Nt = params.Nt; Nr = params.Nr;
    Lt = params.Lt; Lr = params.Lr;
    K = params.K; M = size(F_train, 1);
    Gt = params.Gt; Gr = params.Gr;

    % -- 1. Construct the Dictionary Matrix Psi
    % Create virtual transmit and receive steering vectors for the grid
    At_dict = zeros(Nt, Gt);
    Ar_dict = zeros(Nr, Gr);
    angle_grid = linspace(0, pi, Gt); % Assuming Gt = Gr

    for g = 1:Gt
        At_dict(:, g) = (1/sqrt(Nt)) * exp(1j * pi * (0:Nt-1)' * cos(angle_grid(g)));
    end
    for g = 1:Gr
        Ar_dict(:, g) = (1/sqrt(Nr)) * exp(1j * pi * (0:Nr-1)' * cos(angle_grid(g)));
    end
    % The dictionary Psi is the Kronecker product of the conjugated virtual
    % transmit dictionary and the virtual receive dictionary [R-4, eq. (11)]
    Psi = kron(conj(At_dict), Ar_dict);

    % -- 2. Construct the Sensing Matrix Phi
    Phi = zeros(M * Lr, Nt * Nr);
    for m = 1:M
        Fm = F_train{m};
        Wm = W_train{m};
        qm = q_train{m};
        row_start = (m-1)*Lr + 1;
        row_end = m*Lr;
        
        % The Kronecker product was causing a dimension mismatch error.
        % This has been replaced with a manual, explicit construction of the
        % sensing matrix block, which is more robust.
        P_m_T = (Fm * qm).'; % This is a 1 x Nt row vector
        W_m_H = Wm';         % This is an Lr x Nr matrix

        % Manually construct the block of Phi corresponding to the m-th measurement
        Phi_m_block = zeros(Lr, Nt * Nr);
        for i = 1:Nt
            col_start = (i-1)*Nr + 1;
            col_end = i*Nr;
            % Scale the combiner matrix by the i-th element of the precoder vector
            Phi_m_block(:, col_start:col_end) = P_m_T(i) * W_m_H;
        end
        Phi(row_start:row_end, :) = Phi_m_block;
    end

    % -- 3. Generate Received Measurement Vectors and True Sparse Vectors
    y_vec = cell(K, 1);
    h_v = cell(K, 1);

    for k = 1:K
        Hk = H_freq{k};
        y_k_stacked = zeros(M * Lr, 1);
        
        % Generate combined noise for the entire stacked vector
        noise_block = sqrt(noise_var/2) * (randn(M*Lr, 1) + 1j*randn(M*Lr, 1));

        for m = 1:M
            Fm = F_train{m};
            Wm = W_train{m};
            qm = q_train{m};
            tmk = t_train(k, m);

            % --- FIX ---
            % Broke down the matrix multiplication into explicit steps to resolve
            % a dimension mismatch error.
            effective_precoder = Fm * qm; % Should result in an Nt x 1 vector
            channel_times_precoder = Hk * effective_precoder; % Should result in an Nr x 1 vector
            rx_signal_m_no_pilot = Wm' * channel_times_precoder; % Should result in an Lr x 1 vector
            rx_signal_m = rx_signal_m_no_pilot * tmk; % Scale by scalar pilot
            % --- END FIX ---
            
            % Stack the measurements
            row_start = (m-1)*Lr + 1;
            row_end = m*Lr;
            y_k_stacked(row_start:row_end) = rx_signal_m;
        end
        
        % Add combined noise and compensate for pilot symbols
        y_k_stacked = y_k_stacked + noise_block;
        y_k_stacked = y_k_stacked ./ repelem(t_train(k,:)', Lr);

        y_vec{k} = y_k_stacked;

        % Calculate the true sparse vector h_v[k] for this subcarrier
        % h_v[k] = vec(Delta_v[k]) = Psi' * vec(H[k])
        % This projects the true channel onto the dictionary basis
        h_v{k} = Psi' * reshape(Hk, Nt * Nr, 1);
    end
end
