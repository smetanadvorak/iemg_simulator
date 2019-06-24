function [ mfr ] = generate_mu_maxfr( rt, type, par1, par2 )
%[ mfr ] = generate_mu_minimum_fr( rt, type, par1, par2 )
% This function takes recruitment thresholds rt and generates minimum
% (maximum) firing rates for them according to following models:
%1) Fuglevand et al. 1993: same MFR across all the MUs, equal to 8 p/s or
%set by the user in par1; type = 'constant'.
%2) Linear RT model, MFR = par2 + par1 * rt; 
% type = 'linear_rt';
%3) Linear index model, MFR = par2 + par1 * (0:length(rt)-1); type = 'linear_idx';
%4) Keenan and Valero-Cuevas 2007. Random model, MFR is linearly mapped 
% between two random values from U(6,12); so that first value is MFR of 
% smallest MU and second one is MFR of the largest MU. type = 'random';

switch type
    case 'constant'
        if nargin < 3
            par1 = 25; 
        end
        mfr = ones(size(rt)) * par1;
        
    case 'linear_rt'
        mfr = par2 + par1 * rt;
        
    case 'linear_idx'
        mfr = par2 + par1 * (0:length(rt)-1);
        
    case 'random'
        a = 20 + randi(10); %25 + 25 according to the paper (Keenan 2007)
        b = 20 + randi(10); %25 + 25 
        mfr = a + (b-a)/(length(rt)-1) * (0:length(rt)-1);
end

