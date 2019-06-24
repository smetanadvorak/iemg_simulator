function assignment = assign_mf2mn(mf_centers, mu_centers, innerv_area, n_neighbours, muscle_radius)
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

sigma = @(ia) eye(2) * ia;

for mu = 1:nMU
    probfun = @(x,y) reshape(mvnpdf([x(:), y(:)], mu_centers(mu,:), sigma(innerv_area(mu))), size(x)); % adapted to integral2
    in_circle_int = integral2(probfun, -muscle_radius, muscle_radius, borderfun_neg, borderfun_pos);
    out_circle_coeff(mu) = 1/in_circle_int;
end



%% Distmat for supressing closest neighbours
neighbours = knnsearch(mf_centers, mf_centers, 'K', n_neighbours+1);
neighbours = neighbours(:,2:end);



%% Assignment procedure
nMF = size(mf_centers,1);
assignment = nan(nMF,1);
randomized_mf = randperm(nMF);

i = 0;
for mf = randomized_mf
    
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
        clust_hood = mvnpdf(mf_centers(mf,:)', mu_centers(mu,:)', sigma(innerv_area(mu)));
        clust_hood = clust_hood * out_circle_coeff(mu);
        
        % Final a posteriori probability
        probs(mu) = apriori_prob * clust_hood;
    end
    
    probs = probs/sum(probs);
    assignment(mf) = randsample(nMU, 1, true, probs);

    i = i+1;
    if ~mod(i,1000)
        fprintf('%d of muscle fibers assigned\n', i);
    end
end


end

