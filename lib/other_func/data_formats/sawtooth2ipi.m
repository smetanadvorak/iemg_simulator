function [ ipi_filled, ipi_filled_seamless ] = sawtooth2ipi( sawtooth, ipi_saturation )
    
    if nargin < 2
        ipi_saturation = inf;
    end
    
    spikes = sawtooth2spikes(sawtooth);
    i_sawtooth = spikes2sawtooth(spikes(end:-1:1, :), sawtooth(1,1));
    i_sawtooth = i_sawtooth(end:-1:1, :);
    ipi_filled = sawtooth + i_sawtooth;
    ipi_filled = min(ipi_filled, ipi_saturation);
    
    ipi_filled_seamless = ipi_filled;
    ipi_filled_seamless(ipi_filled == 0) = ipi_filled(find(ipi_filled == 0) - 1); % Sew the zeros with preceding values
    ipi_filled_seamless = min(ipi_filled_seamless, ipi_saturation);
end

