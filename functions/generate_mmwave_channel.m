function [H_freq, At, Ar, alpha, tau, theta_tx, phi_rx] = generate_mmwave_channel(params)
    % generate_mmwave_channel - Creates a frequency-selective mmWave MIMO channel.
    %
    % Syntax: [H_freq, At, Ar, alpha, tau, theta_tx, phi_rx] = generate_mmwave_channel(params)
    %
    % Inputs:
    %   params - A struct containing simulation parameters from initialize_parameters().
    %
    % Outputs:
    %   H_freq   - K x 1 cell array, where each cell contains the Nr x Nt channel matrix for a subcarrier.
    %   At       - Nt x L matrix of transmit array steering vectors.
    %   Ar       - Nr x L matrix of receive array steering vectors.
    %   alpha    - L x 1 vector of complex path gains.
    %   tau      - L x 1 vector of path delays.
    %   theta_tx - L x 1 vector of angles of departure (AoD).
    %   phi_rx   - L x 1 vector of angles of arrival (AoA).
    %
    % Description:
    %   This function implements the geometric channel model from [R-4, eq. (2)].
    %   It generates random path delays, angles, and gains to construct the
    %   time-domain and frequency-domain channel matrices.

    % Extract parameters for convenience
    Nt = params.Nt;
    Nr = params.Nr;
    L = params.L;
    Nc = params.Nc;
    K = params.K;
    Ts = params.Ts;
    rolloff = params.rolloff;

    % -- 1. Generate Random Channel Parameters
    % Path gains: Complex Gaussian distribution
    alpha = (randn(L, 1) + 1j * randn(L, 1)) / sqrt(2);

    % Path delays: Uniformly distributed
    tau = rand(L, 1) * (Nc - 1) * Ts;

    % Angles of Departure (AoD) and Arrival (AoA): Uniformly distributed in [0, 2*pi]
    % For ULA, cos(angle) is uniform in [-1, 1], which means angle is not uniform.
    % We will assume angles are uniform in [0, pi] for simplicity as per common practice.
    theta_tx = pi * rand(L, 1); % AoD
    phi_rx   = pi * rand(L, 1); % AoA

    % -- 2. Create Array Steering Vectors
    At = zeros(Nt, L);
    Ar = zeros(Nr, L);
    for l = 1:L
        At(:, l) = (1/sqrt(Nt)) * exp(1j * pi * (0:Nt-1)' * cos(theta_tx(l)));
        Ar(:, l) = (1/sqrt(Nr)) * exp(1j * pi * (0:Nr-1)' * cos(phi_rx(l)));
    end

    % -- 3. Create Time-Domain Channel Taps (H_d)
    H_time = zeros(Nr, Nt, Nc);
    for d = 0:(Nc-1)
        H_d = zeros(Nr, Nt);
        for l = 1:L
            % Raised-cosine pulse shaping filter effect
            t = d * Ts - tau(l);
            if t == 0
                p_rc = 1;
            elseif abs(t) == Ts / (2 * rolloff)
                p_rc = (pi/4) * sinc(1/(2*rolloff));
            else
                p_rc = sinc(t/Ts) * cos(pi * rolloff * t / Ts) / (1 - (2 * rolloff * t / Ts)^2);
            end
            H_d = H_d + alpha(l) * p_rc * (Ar(:,l) * At(:,l)');
        end
        % Normalization factor from [R-4, eq. (2)]
        H_time(:,:,d+1) = sqrt(Nt * Nr / L) * H_d;
    end

    % -- 4. Convert to Frequency-Domain Channel (H[k])
    H_freq = cell(K, 1);
    for k = 0:(K-1)
        H_k = zeros(Nr, Nt);
        for d = 0:(Nc-1)
            H_k = H_k + H_time(:,:,d+1) * exp(-1j * 2 * pi * k * d / K);
        end
        H_freq{k+1} = H_k;
    end
end
