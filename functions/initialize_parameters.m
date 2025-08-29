function params = initialize_parameters()
    % initialize_parameters - Sets up simulation parameters for the mmWave project.
    %
    % Syntax: params = initialize_parameters()
    %
    % Outputs:
    %   params - A struct containing all the simulation parameters.
    %
    % Description:
    %   This function defines all the necessary parameters for the mmWave
    %   channel estimation project, based on the specifications in [R-4].
    %   Using a single function to manage parameters makes the main simulation
    %   scripts cleaner and easier to modify.

    % -- System Parameters
    params.Nt = 32;             % Number of transmit antennas
    params.Nr = 32;             % Number of receive antennas
    params.Lt = 1;              % Number of transmit RF chains
    params.Lr = 4;              % Number of receive RF chains
    params.Ns = 1;              % Number of data streams (equal to Lt)
    params.K = 16;              % Number of OFDM subcarriers

    % -- Channel Parameters
    params.L = 4;               % Number of channel paths
    params.Nc = 4;              % Number of channel delay taps
    params.Ts = 1/1760e6;       % Sampling period (from IEEE 802.11ad)
    params.rolloff = 0.8;       % Rolloff factor for the raised-cosine filter

    % -- Compressive Sensing Parameters
    params.Gt = 64;             % Size of the transmit dictionary (quantized angles)
    params.Gr = 64;             % Size of the receive dictionary (quantized angles)
    params.Nq = 2;              % Number of quantization bits for phase shifters

    % -- Simulation Control
    % SNR range for NMSE vs. SNR plot
    params.SNR_dB_range_nmse = -15:5:10;
    % Training frames for NMSE vs. SNR plot
    params.M_values_nmse = [80, 120];

    % Training frames range for NMSE vs. M plot
    params.M_range_frames = 20:20:100;
    % SNR values for NMSE vs. M plot
    params.SNR_dB_values_frames = [-10, -5, 0];

    % SNR range for Spectral Efficiency plot
    params.SNR_dB_range_spec = -15:5:10;
    % Training frames for Spectral Efficiency plot
    params.M_spec = 60;

    % Number of Monte Carlo trials for averaging results
    params.num_trials = 100; % A smaller number for faster testing, increase for final results

end
