% seed_rng.m
% Set RNG seeds for reproducibility across MATLAB versions.
function seed_rng(seed)
if nargin<1 || isempty(seed), seed = 233; end
rng(seed,'twister');
end
