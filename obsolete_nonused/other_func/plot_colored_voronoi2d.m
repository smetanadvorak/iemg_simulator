function [ ax ] = plot_colored_voronoi2d(pts, labels, Rmuscle)

[v,c] = voronoin(pts, {'Qbb'});
n_classes = max(labels);
colors = distinguishable_colors(n_classes); % https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors  

figure; ax = axes;
hold on;
for cl = 1:n_classes
    inds = find(labels==cl);
    for i = 1:numel(inds)
        vertices = [v(c{inds(i)},1),v(c{inds(i)},2)];
        % Project vertices outside the muscle region onto the muscle border
        vertices_norms = sqrt(vertices(:,1).^2 + vertices(:,2).^2);
        for j = 1:numel(vertices_norms)
            if vertices_norms(j) > Rmuscle
                vertices(j,:) = vertices(j,:) * Rmuscle / vertices_norms(j);
            end
        end
        fill(vertices(:,1), vertices(:,2), colors(cl,:));
    end
end
hold off;

xlim([-Rmuscle, Rmuscle]); ylim([-Rmuscle, Rmuscle]);
axis equal;  axis([-Rmuscle, Rmuscle, -Rmuscle, Rmuscle]);
xlabel('x, mm'); ylabel('y, mm');

end

