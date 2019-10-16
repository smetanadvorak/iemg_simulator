function [ diam_means, diam_stds ] = get_mf_dist_parameters( innervation_areas )

std_d = 9e-6;   % 9 um, value taken from HAMILTON-WRIGHT AND STASHUK: 
                    %PHYSIOLOGICALLY BASED SIMULATION OF CLINICAL EMG SIGNALS   
%std_d = 6e-6;    % what I get from Voronoi                   
mean_d = 55e-6; % 55 um, taken from the same source
cv = std_d/mean_d;

innervation_areas = innervation_areas/sum(innervation_areas);
cumul_ia = cumsum(innervation_areas);
Nmu = numel(innervation_areas);
diam_means = zeros(Nmu, 1);
diam_stds = zeros(Nmu, 1);

%% A method by Stashuk and Hamilton-Wright (PHYSIOLOGICALLY BASED SIMULATION OF CLINICAL EMG SIGNALS)
for i = 1:Nmu
%    diam_means(i) = mean_d - std_d + innervation_areas(i) * 2*std_d; % See eq. (9) in the article
    diam_means(i) = mean_d - std_d + cumul_ia(i) * 2*std_d; % I don't know, but the original function worked bad for large number of MUs
    diam_stds(i) = diam_means(i) * cv; % See eq. 11
end

end

