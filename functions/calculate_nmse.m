function nmse = calculate_nmse(H_true_freq, H_est_freq)
    % calculate_nmse - Computes the Normalized Mean Square Error between two channels.
    %
    % Syntax: nmse = calculate_nmse(H_true_freq, H_est_freq)
    %
    % Inputs:
    %   H_true_freq - K x 1 cell array of the true frequency-domain channel matrices.
    %   H_est_freq  - K x 1 cell array of the estimated frequency-domain channel matrices.
    %
    % Outputs:
    %   nmse - The calculated NMSE value in dB.
    %
    % Description:
    %   This function implements the NMSE formula from [R-4, eq. (34)]. It sums
    %   the squared error and the channel energy across all subcarriers.

    K = length(H_true_freq);
    total_error_power = 0;
    total_channel_power = 0;

    for k = 1:K
        H_true = H_true_freq{k};
        H_est  = H_est_freq{k};

        % Sum of squared Frobenius norms for error
        total_error_power = total_error_power + norm(H_est - H_true, 'fro')^2;

        % Sum of squared Frobenius norms for channel energy
        total_channel_power = total_channel_power + norm(H_true, 'fro')^2;
    end

    % Calculate NMSE
    nmse_linear = total_error_power / total_channel_power;
    
    % Convert to dB
    nmse = 10 * log10(nmse_linear);
end
