function spikes = sawtooth2spikes( sawtooth )
% spikes = sawtooth2spikes(sawtooth)

nMU = min(size(sawtooth));

if size(sawtooth,1) < size(sawtooth,2)
    sawtooth = sawtooth';
end

spikes = zeros(max(size(sawtooth)), nMU);

for i = 1:nMU
    spikes(sawtooth(:,i) == 0, i) = 1;
end

