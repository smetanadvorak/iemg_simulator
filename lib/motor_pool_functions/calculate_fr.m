function [ fr ] = calculate_fr(rt, minfr, maxfr, frs, excitation)

if excitation > rt
    fr = minfr + max(0, frs*(excitation - rt));
    fr = min(maxfr, fr);
else
    fr = 0;
end

end

