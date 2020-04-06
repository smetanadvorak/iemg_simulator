function [shifted] = hr_shift_template(x, delay)
% This function takes waveform 'x' and shifts it by a subsample step
% 'delay'. This is based on assumption that fft of x is the same as fft of
% the real form. This function linearly transforms the phase of fft(x) and then
% gives back ifft(X_shifted);
% Delay is measured in sampling period, so that delay = 0.1 means 1/10th of the
% sampling period.

x = x(:);

if ~mod(numel(x),2)
    x = [x; 0];
    padded = 1;
else
    padded = 0;
end

N = length(x);

X = fft(x);
X0 = X(1);
Xk = X(2:ceil(end/2));

Sk = Xk .* exp(1i * (2*pi*delay) * (1:length(Xk))'/N);
S = [X0; Sk; conj(Sk(end:-1:1))];
shifted = ifft(S);

if padded
    shifted = shifted(1:end-1);
end