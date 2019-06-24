function Vm = get_tm_current(z, D1, D2)
% This transmembrane current model is taken from:
% P. Rosenfalck Intra and extracellular fields of active nerve and muscle fibers. 
% A physico-mathematical analysis of different models, 1969;
% And from:
% D. Farina and R. Merletti. A novel approach for precise simulation of 
% the EMG signal detected by surface electrodes, 2001;
% Vm = get_sfap(z, D1, D2); z is in millimeters, Vm is in millivolts.

if nargin < 3
    D2 = -90; % -90 mV
end
if nargin < 2
    D1 = 96; % 96 mv/mm^3
end

Vm = nan(size(z));
Vm(z >  0) = D1 * (z(z > 0).^3) .* exp(-z(z > 0)) + D2;
Vm(z <= 0) = D2;

end

