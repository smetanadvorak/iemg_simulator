function [ spikes ] = firings2spikes( firings, len, fs )
nMU = length(firings);
spikes = zeros(len, nMU);

if nargin < 3
    fs = 1;
end

firings_cropped = 0;

for i=1:nMU
    if ~isempty(firings{i}) && firings{i}(end)*fs > len
        firings_cropped = 1;
        firings{i}(firings{i}*fs > len) = [];
    end
    spikes(round(firings{i}*fs), i) = 1;
end

if firings_cropped
    warning('firings2spikes: firings exceed the specified length, result is cropped to fit. Not a big deal, but check if needed.');
end
end

