function output = spikes2firings(spikes, fs)
%firings = spikes2firings( spikes, fs ) returns firings (nMU
%cells with firing times) from given spikes (T x nMU vector of zeros and ones
%where ones designate a firing)
%Also accepts 3-dimensional spikes matrix (3rd dimension for channel) if
%different spike trains are available for different channels. In that case
%output is a cell(n_channels, 1) of cells(nMU,1);
%
if nargin < 2
    fs = 1;
end

len = size(spikes, 1);
nMU = size(spikes, 2);
nCH = size(spikes, 3); % For Multichannel spikes

firings = cell(nMU,1);
firings_multichan = cell(nCH,1);

for ch = 1:nCH
    for m = 1:nMU
        tmp = spikes(:, m, ch); % Initial bad indexing order -> speed up
        for t = 1:len
            if tmp(t)
                firings{m} = [firings{m}; t/fs];
            end
        end
    end
    firings_multichan{ch} = firings;
end
    
if nCH > 1
    output = firings_multichan;
else
    output = firings;
end

end

