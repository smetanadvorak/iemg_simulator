function [ profile ] = generate_identification_input( len, pause, levels, amplitude )
% profile = generate_identification_input( len, pause, levels, amplitude )
% This function generates the identification profiles for the muscle force
% prediction experiments. It approximates three main classes of the dayly
% contraction profiles: steep contractions until desired level, moderately
% increasing/decreasing trapezoidal contractions and slow variations around
% a set level (contact force adjustment).
%
% len sets the lengths of each profile, if a vector of structure [rect,
% trapez, const], defines different lengts for each profile. If scalar x,
% sets the lengths to rect = x, trapez = 5x, const = 5x.
%
% pause sets the pause between profiles. Just ensure that it's longer than
% full relaxation time of the muslce (around 250 ms should be sufficient)
%
% levels is a vector with the plateau excitation levels;
% 
% amplitude sets the variation amplitude for the constant profile, relative
% to the set level. Default is 0.2.

if nargin < 4
    amplitude = 0.1;
end

if numel(len)==1
    rect_len = len;
    trapez_len = 5*len;
    const_len = 4*len;
else
    rect_len = len(1);
    trapez_len = len(2);
    const_len = len(3);
end

profile = [];

% Part 1: rectangles
for i = 1:numel(levels)
    profile = [profile; levels(i) * ones(rect_len,1); zeros(pause, 1)];
end
    
% Part 2: trapezoids
profile_nodes = [0, 40, 60, 100]/100 * trapez_len; % trapezoidal profile
profile_timeline = linspace(0,trapez_len, trapez_len)';
for i = 1:numel(levels)
    profile_vals = [0, levels(i) + eps, levels(i) + eps, 0]; % trapezoidal profile
    pr = interp1(profile_nodes, profile_vals, profile_timeline(:));
    profile = [profile; pr; zeros(pause,1)];
end

%Part 3: constant with slow oscillations
profile_nodes = [0, 20, 100]/100 * const_len; 
profile_timeline = linspace(0,const_len, const_len)';
for i = 1:numel(levels)
    profile_vals = [0, levels(i) + eps, levels(i) + eps]; 
    pr = interp1(profile_nodes, profile_vals, profile_timeline);
    pr = pr + amplitude * levels(i) * cos(2*pi/(const_len/4) * linspace(0,const_len, const_len)') - amplitude * levels(i);
    profile = [profile; pr; zeros(pause,1)];
end

end

