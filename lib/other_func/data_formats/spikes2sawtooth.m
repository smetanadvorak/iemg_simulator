function [ seq ] = spikes2sawtooth(spikes, initial_values)

if nargin < 2
    initial_values = ones(1, size(spikes,2));
end

l = size(spikes, 1);
w = size(spikes, 2);

seq = zeros(l,w);

initial_values = initial_values .* (spikes(1,:) ~= 1);
seq(1,:) = initial_values;

for j = 1:w
    for i = 2:l
        if spikes(i,j)
            seq(i,j) = 0;
        else
            seq(i,j) = seq(i-1,j) + 1;
        end
    end
end

end

