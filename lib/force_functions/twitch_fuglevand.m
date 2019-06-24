function [ twitch ] = twitch_fuglevand(t, P, T)
% This function implements the twitch model presented in Fuglevand et al. -
% Models of Recruitment and Rate coding organisation in motor units, 1993.

twitch = P/T * t .* exp(1 - t/T); 


end

