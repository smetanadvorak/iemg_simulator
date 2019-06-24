

analysisN = 99;

[idx, cluster_centers] = kmeans(mu_pool.mf_centers(mu_pool.assignment == analysisN, :), n_branches(analysisN));

figure; 
axis tight;

plot(mu_pool.mn_pool.centers(analysisN,1), mu_pool.mn_pool.centers(analysisN,2), 'b.', 'markersize', 20); hold on;
xlabel('x, mm'); ylabel('y, mm');

assgnt = find(mu_pool.assignment == analysisN);
for c = 1:size(cluster_centers,1)
    
    plot(mu_pool.mf_centers(assgnt(idx==c), 1), mu_pool.mf_centers(assgnt(idx==c), 2), '.', 'markersize',10);
    plot(cluster_centers(c,1) , cluster_centers(c,2), 'ko');
    line([mu_pool.mn_pool.centers(analysisN,1),cluster_centers(c,1)] , [mu_pool.mn_pool.centers(analysisN,2),cluster_centers(c,2)], 'color', 'k');
    %text(cluster_centers(c,1)+[mu_pool.mn_pool.centers(analysisN,1) - cluster_centers(c,1)]/2 , cluster_centers(c,2)+[mu_pool.mn_pool.centers(analysisN,2) - cluster_centers(c,2)]/2, '$d_k$', 'fontsize', 20);
end

axs = axis;
phi_circle = linspace(0, 2*pi, 1000)';
phi_circle = phi_circle(1:end-1);
muscle_border = [Rmuscle * cos(phi_circle), Rmuscle * sin(phi_circle)];
plot(muscle_border(:,1), muscle_border(:,2), 'k');
axis(axs);

legend('MU innervation center', 'Muscle fibers', 'Branch innervation centers', 'Muscle border');

clear analysisN idx cluster_centers muscle_border c