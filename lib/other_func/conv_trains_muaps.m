function [ emg ] = conv_trains_muaps(trains, muaps, mode)

if nargin < 3
    mode = 'causal';
end

if size(trains,2) ~= size(muaps, 2)
    error('Error in trains or muaps matrices size: expected number of motor units didn''t coincide');
end

emg = 0;

switch mode
    case 'centered'
        for mu = 1:size(trains, 2)
            emg = emg + conv(trains(:,mu), muaps(:,mu), 'same');
        end
    case 'causal'
        for mu = 1:size(trains, 2)
            emg = emg + conv(trains(:,mu), muaps(:,mu), 'full');
        end
        emg = emg(1:end-length(muaps)+1, :);
end

emg(end:size(trains,1)) = 0;

