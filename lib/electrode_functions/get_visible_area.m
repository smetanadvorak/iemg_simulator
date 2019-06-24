function [ visible_area ] = get_visible_area(mf_pts, el_pt, normal)

el_pt = el_pt(:)';
Np = size(mf_pts,1);

% Normal vector of the electrode element (repmat for vectorized calc)
normal_norm = sqrt(sum(normal.^2));
normal = repmat(normal, Np, 1);

% Directions from the electrode element to the mf's elements:
dir_vectors = mf_pts - repmat(el_pt, Np, 1);
dir_norms = sqrt(sum(dir_vectors.^2, 2));
% Normalize them:
% Will be done after taking scalar product ...

% Scalar products of surface's normal and field's radial vectors
% and ormalization by radial vector length. No the
scalar_prods = sum(normal .* dir_vectors, 2) ./ dir_norms / normal_norm ;  

% Normal vectors are not normalized because their length is equal to
% surface area.

scalar_prods = max(0, scalar_prods); % Since electrode is a closed body, we assume 
                                     % no mf-induced potential inside of
                                     % it, which would make (normal, direction) < 0
                                     
% See my report for that square; normal_norm is equal to element's area   
visible_area = (scalar_prods.^2) * normal_norm; 
    
