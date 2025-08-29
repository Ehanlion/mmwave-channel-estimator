function [h_v_est, T_est] = swomp(params, y_vec, Phi, Psi, noise_var)
    % swomp - Implements the Simultaneous Weighted Orthogonal Matching Pursuit algorithm.
    %
    % Syntax: [h_v_est, T_est] = swomp(params, y_vec, Phi, Psi, noise_var)
    %
    % Inputs:
    %   params    - Struct with simulation parameters.
    %   y_vec     - K x 1 cell array of stacked measurement vectors.
    %   Phi       - Sensing matrix.
    %   Psi       - Dictionary matrix.
    %   noise_var - Variance of the AWGN.
    %
    % Outputs:
    %   h_v_est - K x 1 cell array of estimated sparse channel vectors.
    %   T_est   - Vector containing the indices of the estimated common support.
    %
    % Description:
    %   This function implements the SW-OMP algorithm as described in
    %   [R-4, Alg. 1] and the tutorial slides. It iteratively finds the common
    %   sparse support of the channel across all subcarriers.

    % Extract parameters
    K = params.K;
    M = size(y_vec{1}, 1) / params.Lr;
    Lr = params.Lr;
    max_paths = params.L * 2; % Stop after finding more paths than expected

    % -- Step 0: Initialization
    % Combine sensing and dictionary matrices
    Upsilon = Phi * Psi;
    [~, N_total] = size(Upsilon);

    % Initialize residuals to the input signals
    r = y_vec;
    T_est = []; % Estimated support
    h_v_est = cell(K, 1);
    for k=1:K
        h_v_est{k} = zeros(N_total, 1);
    end

    % Stopping criterion: residual energy vs. noise energy
    % The total noise power in the stacked vector is M * Lr * noise_var
    stop_threshold = M * Lr * noise_var;
    
    current_mse = inf;
    iter = 0;

    % -- Main Iteration Loop
    while current_mse > stop_threshold && iter < max_paths
        iter = iter + 1;

        % -- Step 1: Correlation
        correlations = zeros(N_total, 1);
        for k = 1:K
            % Project residual onto the dictionary atoms
            c_k = Upsilon' * r{k};
            % Sum the magnitudes across all subcarriers
            correlations = correlations + abs(c_k);
        end

        % -- Step 2: Find the best candidate atom
        [~, p_star] = max(correlations);
        
        % Update the common support
        T_est = union(T_est, p_star);

        % -- Step 3: Project measurements onto the new support (LS solution)
        Upsilon_T = Upsilon(:, T_est);
        for k = 1:K
            % Solve the least-squares problem: x = (A'A)^-1 * A' * y
            x_k_T = pinv(Upsilon_T) * y_vec{k};
            
            % Place the estimated gains into the full sparse vector
            h_v_est{k}(T_est) = x_k_T;
        end

        % -- Step 4: Update the residual
        total_residual_power = 0;
        for k = 1:K
            r{k} = y_vec{k} - Upsilon * h_v_est{k};
            total_residual_power = total_residual_power + norm(r{k})^2;
        end
        
        % -- Step 5: Check stopping criterion
        current_mse = total_residual_power / K;
    end
end
