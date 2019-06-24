function [assignment, adopted_by_mu] = assign_mf_gauss_size(mf_centers, mu_centers, innerv_area, diam_means, diam_stds, mf_diameters, muscle_radius, intermingling_factor)
% This function takes a set of motor neurons with their innervation centers
% and areas, as well as coordinates of muscle fibers and assigns fibers to
% neurons according to several rules:
% Proximity to the innervation center;
% Innervation area (or target number of innervated fibers)
% Self-avoiding phenomena (see n_neighbours variable)


%% Out-of-muscle area compensation
% Calculates how much of the MU's gaussian distribution is outside of the
% muscle border and inflates the rest of the distribution according to it
nMU = size(mu_centers,1);
borderfun_pos = @(x)real( sqrt(muscle_radius^2 - x.^2));
borderfun_neg = @(x)real(-sqrt(muscle_radius^2 - x.^2));
out_circle_coeff = ones(nMU,1);

for mu = 1:nMU
    probfun = @(x,y) reshape(mvnpdf([x(:), y(:)], mu_centers(mu,:), eye(2)*innerv_area(mu)), size(x)); % adapted to integral2
    in_circle_int = integral2(probfun, -muscle_radius, muscle_radius, borderfun_neg, borderfun_pos);
    out_circle_coeff(mu) = 1/in_circle_int;
end



%% Intermingling factor adjusts the MUs' areas by inflating their covariation matrices
% I no longer use this, set coefficient to 1
intermingling_coeff = 1; %1/chi2inv(1 - intermingling_factor, 2);
% sigma = percentile_radius/epsilon = mu_area/intermingling_coeff
% http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
inflated_std = innerv_area * intermingling_coeff;



%% Distmat for supressing closest neighbours
n_neighbours = 3;
neighbours = knnsearch(mf_centers, mf_centers, 'K', n_neighbours+1);
neighbours = neighbours(:,2:end);



%% Assignment procedure

nMF = size(mf_centers,1);

expected_per_mu = ceil(innerv_area/sum(innerv_area) * nMF);
adopted_by_mu = zeros(size(expected_per_mu));
assignment = nan(nMF,1);

randomized_mf = randperm(nMF);
i = 0;
for mf = randomized_mf
    %sum_percent_adopted = sum((expected_per_mu - adopted_by_mu)./expected_per_mu);
    %sum_percent_adopted = sum_percent_adopted/nMU;
    probs = zeros(nMU,1);
    
    for mu = 1:nMU
        % Supression assignment if neighbours are from the same MU
        % Promotes intermingling
        % Problem: spreads the MUs too much;
        if any(assignment(neighbours(mf,:)) == mu)
            probs(mu) = 0;
            continue;
        end

        % A priori probability of the assignment
        apriori_prob = innerv_area(mu)/sum(innerv_area);
        
        % Likelihood coming from clustered nature of mf distribution
        clust_hood = mvnpdf(mf_centers(mf,:)', mu_centers(mu,:)', eye(2)*inflated_std(mu));
        clust_hood = clust_hood * out_circle_coeff(mu);
        %clust_hood = clust_hood / mvnpdf([0;0], [0;0], eye(2)*inflated_std(mu));
        
        % Likelihood coming from the average mf radius of this MU
        %rad_hood = normpdf(mf_diameters(mf), diam_means(mu), diam_stds(mu));
        rad_hood = 1;    
        
        % Final a posteriori probability
        probs(mu) = apriori_prob * clust_hood * rad_hood;
    end
    
    assignment(mf) = randsample(nMU, 1, true, probs/sum(probs));
    adopted_by_mu(assignment(mf)) = adopted_by_mu(assignment(mf)) + 1;

    i = i+1;
    if ~mod(i,1000)
        fprintf('%d of muscle fibers assigned\n', i);
    end
end


end

