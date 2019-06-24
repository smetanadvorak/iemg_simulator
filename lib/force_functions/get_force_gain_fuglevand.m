function gain = get_force_gain_fuglevand( ipi, T )
% gain = get_force_gain_fuglevand( ipi, T ) returns the gain value for the
% force output for a motor unit with current inter-pulse-interval ipi and
% T-parameter of the twitch T. This function corresponds to Fuglevands
% nonlinear gain model for the force output, see Fuglevant - Models of Rate
% Coding ..., eq. 17.

Sf = @(x)(1 - exp( -2*(x).^3));

inst_dr = T./ipi;
gain = ones(size(inst_dr));
gain(inst_dr > 0.4) = Sf(inst_dr(inst_dr > 0.4)) ./ inst_dr(inst_dr > 0.4) / (Sf(0.4)/0.4);  % See Fuglevand, eq. 17

