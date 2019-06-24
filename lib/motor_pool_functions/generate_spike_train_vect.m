function [ spikes, next_state, last_ipi] = generate_spike_train_vect( t, prev_state, excitation, rt, minfr, maxfr, frs, CV, fs )
    T = numel(t);
    N = numel(rt);
    spikes = zeros(T,N);
    ipi = zeros(N);
    next_firing = prev_state;
    next_state = zeros(N,1);
    last_ipi = zeros(N,1);
    
    for n = 1:N
        for i = 1:T
            if excitation(i) > rt(n)
                if isnan(next_firing(n))
                    ipi(n) = fs/calculate_fr(rt(n), minfr(n), maxfr(n), frs(n), excitation(i)); % mean IPI
                    ipi(n) = ipi(n) + randn * ipi(n) * CV; % add some variation
                    next_firing(n) = t(i) + round(ipi(n));
                end
                if t(i) == next_firing(n)
                    spikes(i,n) = 1;
                    ipi(n) = fs/calculate_fr(rt(n), minfr(n), maxfr(n), frs(n), excitation(i));
                    ipi(n) = ipi(n) + randn * ipi(n) * CV;
                    next_firing(n) = t(i) + round(ipi(n));
                end
            else
                next_firing(n) = nan;
            end
        end
        last_ipi(n) = ipi(n);
        next_state(n) = next_firing(n);
    end
end

