function [diameters] = assign_mf_diameters(mf_assignment, diam_means, diam_stds)

diameters = zeros(size(mf_assignment));
Nmu = max(mf_assignment);

%% A method by Stashuk and Hamilton-Wright (PHYSIOLOGICALLY BASED SIMULATION OF CLINICAL EMG SIGNALS)

for i = 1:Nmu
    ind = find(mf_assignment==i);
    diameters(ind) = diam_means(i) + diam_stds(i) * randn(size(ind));
end
