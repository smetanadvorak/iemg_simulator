function [ frs ] = generate_mu_fr_slope( rt, type, par1, par2 )
%[ mfr ] = generate_mu_fr_slope( rt, type, par1, par2 )
% This function takes recruitment thresholds rt and generates firing rate
% slopes to following models:
%1) Fuglevand et al. 1993: same FRS across all the MUs, equal to par1
%[p/(s*excitation)];
% ; type = 'constant'.
%2) Linear RT model, FRS = par2 + par1 * rt; 
% type = 'linear_rt';
%3) Linear index model, FRS = par2 + par1 * (0:length(rt)-1); type = 'linear_idx';
%4) Random model, FRS is linearly mapped between two random values from U(6,12); 
% so that first value is MFR of smallest MU and second one is MFR of the largest MU. type = 'random';
% This model is analogous to that of Keenan and Valero-Cuevas 2007 for
% maxfr.

switch type
    case 'constant'
        if nargin < 3
            par1 = 8; 
        end
        frs = zeros(size(rt)) + par1;
        
    case 'linear_rt'
        frs = par2 + par1 * rt;
        
    case 'linear_idx'
        frs = par2 + par1 * (0:length(rt)-1);
        
    case 'random'
        a = 5 + randi(7);
        b = 5 + randi(7); 
        frs = a + (b-a)/(length(rt)-1) * (0:length(rt)-1);
end

