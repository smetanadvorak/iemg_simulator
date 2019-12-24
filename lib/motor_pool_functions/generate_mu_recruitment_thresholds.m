function [rt,rtz] = generate_mu_recruitment_thresholds( N, RR, RM, b, mode )
%function [rt,rtz] = generate_mu_recruitment_thresholds( N, RR, RM, b, mode )
%
%Generates N motor neurons' threshold values. This model is a slight variation
%of the one proposed in Fuglevand 1993 (Models of Recruitment and Rate Coding ...).
%Their formula is: rt(i) = exp(i * ln(RR) / N).
%
%The idea was developed in De Luca (Hierarchical control of motor units
%...). Their formula is rt(i) = 1/100 * exp(ln(RR)/N * i).
%
%My formula adjusts both the recruitment range (RR = rt(N)/rt(1)) and the
%value of maximum RM = rt(N).
%
%N is number of motor neurons;
%RR is range of recruitment thresholds (rt(N)/rt(1)).
%RM is threshold of the largest motor neuron;
%mode is either 'akhmadeev', 'fuglevand' or 'deluca'.
%b is slope correction coefficient for De Luca's model
%
% mode:
%'fuglevand': uses only N and RR;
%'deluca': uses only N, RR and b;
%'ls2n': uses N, RR, RM;
%
%rt is recruitment thresholds
%rtz is almost the same, but starts from zero, convenient in simulation;

if nargin < 5
    mode = 'ls2n';
end

i = 1:N;

switch mode
    case 'fuglevand'
        rt = exp(i*ln(RR)/N);
        rtz = rt - rt(1);
        rtz = rtz * max(rt)/max(rtz);
        
    case 'deluca'
        rt = (b*i/N)/100 .* exp(i*ln(RR/b));
        rtz = rt - rt(1);
        rtz = rtz * max(rt)/max(rtz);
        
    case 'ls2n'
        rt = (RM/RR) * exp((i-1) * log(RR) / (N-1));
        % Version starting from zero
        rtz = (RM/RR) * (exp((i-1) * log(RR+1) / N) - 1);
        rtz = rtz * max(rt)/max(rtz);
        
end

end

