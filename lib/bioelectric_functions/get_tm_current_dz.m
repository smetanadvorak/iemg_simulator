function Vm = get_tm_current_dz(z, D1)
% This first derivative of transmembrane current model taken from:
% P. Rosenfalck Intra and extracellular fields of active nerve and muscle fibers. 
% A physico-mathematical analysis of different models, 1969;
% And from:
% D. Farina and R. Merletti. A novel approach for precise simulation of 
% the EMG signal detected by surface electrodes, 2001;
% Vm = get_sfap(z, D1); z is in millimeters, Vm is in millivolts.

if nargin < 2
    D1 = 96; % 96 mv/mm^3
end

zp = z(z > 0);
Vm = nan(size(z));
Vm(z >  0) = D1 * (3*(zp.^2) - (zp.^3)).* exp(-zp);
Vm(z <= 0) = 0;

end

