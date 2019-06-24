function vertices = fp_sampling_wrapper( np, density_map, init_p )
% vertices = fp_sampling_wrapper( np, density_map, init_p )
% Perform fast marching random samping of a region described in
% 'density_map'. Points outside of the region should be filled with zeros
% (zero density).
% Wraps around 

if nargin < 3
    init_p = size(density_map)/2;
end

vertices = zeros(2,np);
vertices(:,1) = init_p(:);
for i = 2:np
    D = perform_fast_marching(1./density_map, vertices(:, 1:i-1));
    [~,ind] = max(D(:));
    [x,y] = ind2sub([np np],ind);
    vertices(:,i) = [x;y];
end


end

