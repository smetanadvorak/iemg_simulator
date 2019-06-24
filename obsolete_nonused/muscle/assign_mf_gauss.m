function assignment = assign_mf_gauss(mf_centers, mu_centers, ia, intermingling_factor)

if nargin < 4
    intermingling_factor = 0.05;
end

intermingling_coeff = chi2inv(1 - intermingling_factor, 2);
% sigma = percentile_radius/epsilon = mu_area/intermingling_coeff
% http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/

nMU = size(mu_centers,1);
nMF = size(mf_centers,1);

expected_per_mu = ceil(ia/sum(ia) * nMF);

adopted_by_mu = zeros(size(expected_per_mu));
assignment = nan(nMF,1);

randomized_mf = randperm(nMF);

i = 0;
for mf = randomized_mf
    
    %sum_percent_adopted = sum((expected_per_mu - adopted_by_mu)./expected_per_mu);
    %sum_percent_adopted = sum_percent_adopted/nMU;
    
    probs = zeros(nMU,1);
    
    for mu = 1:nMU
        clust_prob = mvnpdf(mf_centers(mf,:)', mu_centers(mu,:)', eye(2)*ia(mu)/intermingling_coeff);
        clust_prob = clust_prob / mvnpdf([0;0], [0;0], eye(2)*ia(mu)/intermingling_coeff);
        %adopt_prob = (expected_per_mu(mu) - adopted_by_mu(mu))/expected_per_mu(mu)/sum_percent_adopted;
        adopt_prob = (expected_per_mu(mu) - adopted_by_mu(mu))/expected_per_mu(mu);
        probs(mu) = clust_prob;% * adopt_prob;
    end
    
    assignment(mf) = randsample(nMU, 1, true, probs/sum(probs));
%    adopted_by_mu(assignment(mf)) = adopted_by_mu(assignment(mf)) + 1;

    i = i+1;
    if ~mod(i,1000)
        fprintf('%d of muscle fibers assigned\n', i);
    end
end


end

