function [ spikes, next_state, last_ipi] = generate_spike_train( t, prev_state, excitation, rt, minfr, maxfr, frs, CV, fs )
    T = numel(t);
    spikes = zeros(T,1);
    ipi = 0;
    next_firing = prev_state;
    for i = 1:T
        if excitation(i) > rt
            if isnan(next_firing)
                ipi = fs/calculate_fr(rt, minfr, maxfr, frs, excitation(i)); % mean IPI
                ipi = ipi + randn * ipi * CV; % add some variation
                next_firing = t(i) + round(ipi);
            end
            if t(i) == next_firing
                spikes(i) = 1;
                ipi = fs/calculate_fr(rt, minfr, maxfr, frs, excitation(i));
                ipi = ipi + randn * ipi * CV;
                next_firing = t(i) + round(ipi);
            end
        else
            next_firing = nan;
        end
    end
    last_ipi = ipi;
    next_state = next_firing;
end

