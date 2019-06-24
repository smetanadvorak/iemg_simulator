function [ ax ] = plot_colored_voronoi3d(pts, labels, border, ax)

% if nargin < 4
%     ax = axes;
% end

[v,c] = voronoin(pts, {'Qbb'});
n_classes = max(labels);
colors = distinguishable_colors(n_classes); % https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors  


hold on;
for cl = 1:n_classes
    inds = find(labels==cl);
    for i = 1:numel(inds)
        vertices = [v(c{inds(i)},1),v(c{inds(i)},2)];
        z = 2.5 * (cl^1.5) * ones(size(vertices,1),1);
        % Project vertices outside the muscle region onto the muscle border
        vertices_norms = sqrt(vertices(:,1).^2 + vertices(:,2).^2);
        for j = 1:numel(vertices_norms)
            if vertices_norms(j) > border
                vertices(j,:) = vertices(j,:) * border / vertices_norms(j);
            end
        end
        %if ~any(vertices_norms > border)
            fill3(vertices(:,1), vertices(:,2), z, colors(cl,:));
        %end
    end
end
hold off;
view(45,15);
set(gcf, 'position',  [330     1   770   804]);
title('Fibers assignment, layer per MU');

end

