function assignment = assign_mf_hamilton(mf_centers, mu_centers, ia)

nMU = size(mu_centers,1);
nMF = size(mf_centers,1);

expected_per_mu = round(ia/sum(ia) * nMF);
adopted_by_mu = zeros(size(expected_per_mu));
assignment = nan(nMF,1);

weights = [0.6; 0.4; 0.0; 0.0]; % [0.5; 0.3; 0.1; 0.1];
factors = zeros(size(weights));

% Pre-calculate sums of distances from MUs to MFs
%sum_dist = sum(dist(mu_centers, mf_centers'));
sum_dist = mean(dist(mu_centers, mf_centers')); % My modification - normalization by nMU

randomized_mf = randperm(nMF);

i = 0;
for mf = randomized_mf
    % Consider here to choose closest motor units, probably
    
    sum_percent_adopted = sum((expected_per_mu - adopted_by_mu)./expected_per_mu);
    sum_percent_adopted = sum_percent_adopted/nMU; % My modification - normalization
    sum_adopted = sum(adopted_by_mu);
    w = zeros(nMU,1);
    for mu = 1:nMU
%        factors(1) = dist(mf_centers(mf,:), mu_centers(mu,:)') / sum_dist(mf);
        factors(1) = sum_dist(mf) / dist(mf_centers(mf,:), mu_centers(mu,:)');
        factors(2) = (expected_per_mu(mu) - adopted_by_mu(mu))/expected_per_mu(mu)/sum_percent_adopted;
%        factors(3) = adopted_by_mu(mu)/sum_adopted; % Zero at the beginning - weird
%        factors(4) = 1/nMU; % Equal for all MUs - doesn't make sense to me
        factors(3) = 0;
        factors(4) = 0;
        w(mu) = weights' * factors;
    end
    
    [~, assignment(mf)] = max(w);
    
    i = i+1;
    if ~mod(i,1000)
        fprintf('%d of muscle fibers assigned\n', i);
    end
end


end

