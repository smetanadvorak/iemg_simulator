function [ pts, normals, signs, simplified_pts ] = get_needle( len, out_rad, in_rad, circ_segments, lin_segments, tip_segments )

    % Origin is in the center of the tip of the needle
    % Axes are as follows: 
        % - x originates in the center of the tip and is placed
        % along the needle in direction of it's handle.
        % - z looks outwards the needle tip and lays in the plane of it's steepest descent 
        % - y follows right-handed orientation
    % Needle has (all in millimeters)
        % length "len" from the center of the tip to the handle
        % outer radius 
        % inner radius
        % angle of the tip
    % Needle' cannula is modelized as a cylinder, tip is modelised as plane
    
    %% Cylinder approximation
    circle_approx_angles = 2*pi*linspace(0,1,circ_segments+1);
    rad_of_approx_polygon = pi*out_rad / (circ_segments * sin(pi/circ_segments));
    
    circle_approx = [rad_of_approx_polygon * sin(circle_approx_angles(1:end-1));
                     rad_of_approx_polygon * cos(circle_approx_angles(1:end-1))]';
                    % order of sin and cos is inverted because I want the
                    % first point to have y = 0;
    
    dphi = 2*pi / circ_segments;
	dl = len/lin_segments;
    ds_cyl = dl * dphi;
    
    cylinder_points = zeros(lin_segments*circ_segments, 3);
    cylinder_normal = zeros(lin_segments*circ_segments, 3);
    for l = 1:lin_segments
        for c = 1:circ_segments
            cylinder_points((l-1)*circ_segments + c, :) = [dl*(l-1)+dl/2, circle_approx(c,:)];
            cylinder_normal((l-1)*circ_segments + c, :) = [0, circle_approx(c,:)]/norm([0, circle_approx(c,:)]);
            cylinder_normal((l-1)*circ_segments + c, :) = ds_cyl * cylinder_normal((l-1)*circ_segments + c, :);
        end
    end
        
    %% Modelling the tip
    if ~(tip_segments == 5)
        error('This number of tip segments is not implemented yet');
    end
    
    tip_angle = pi/4;
    ellipsis_a = in_rad/sin(tip_angle);
    ellipsis_b = in_rad;
    ds_tip = pi * ellipsis_a * ellipsis_b / 5;
    
    if tip_segments == 5
        center = [0, 0, 0];
        tip_approx_angles = 2*pi*linspace(0,4,5)/4;
        neighbours = [2/3 * ellipsis_a * sin(tip_approx_angles(1:end-1));
                      2/3 * ellipsis_b * cos(tip_approx_angles(1:end-1))]';
                  % order of sin and cos is inverted because I want the
                  % bigger side of ellipsis to be along the z axis
        neighbours = [[cos(tip_angle)*ellipsis_a; 0; -cos(tip_angle)*ellipsis_a; 0], neighbours];
        
        tip = [center; 
               neighbours];
           
        tip_normal = [-cos(tip_angle), 0, sin(tip_angle)];
        tip_normal = tip_normal/norm(tip_normal);
        tip_normal = repmat(tip_normal, tip_segments, 1);
        tip_normal = tip_normal * ds_tip;
        
        simplified_pts = [ neighbours([1,3],:); 
                           neighbours([3,1],:) + repmat([len, 0, 0], 2, 1);
                           neighbours(1,:)];
                       
        
    end
               
    
    pts = [tip; cylinder_points];
    normals = [tip_normal; cylinder_normal];
    signs = [ones(size(tip,1),1); -ones(size(cylinder_points,1),1)];

end

