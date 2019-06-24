classdef IntramuscularArray < Electrode
    
    properties
        step;       % Step of array electrode
    end
    
    methods
        function obj = IntramuscularArray(n_elements, step, differentiation)
            if nargin < 3
                differentiation = 'consecutive';
            end

            obj.n_points = n_elements;
            
            if obj.n_points == 1
                obj.type = 'single point';
                obj.type_short = '1p';
                
            elseif obj.n_points == 2
                obj.type = '2-point differential';
                obj.type_short = '2p_diff';
                
            else
                obj.type = ['intramuscular_', num2str(obj.n_points), 'p_array'];
                obj.type_short = [num2str(obj.n_points), 'p_array'];
            end
                    
            obj.step = step;
            
            % Create the array's points in its own frame (positioned along z-axis)
            obj.pts_origin = [zeros(obj.n_points,2), (0:1:(obj.n_points-1))'*obj.step];
            
            % No normals, since these are point electrodes
            obj.normals_origin = [];%zeros(size(obj.pts_origin,1),3);
            obj.normals_init = [];%zeros(size(obj.pts_origin,1),3);
            obj.normals = [];%zeros(size(obj.pts_origin,1),3);
            
            % Electrode difference matrix
            switch differentiation
                case 'consecutive'
                    obj.diff_mat =          eye(size(obj.pts_origin,1)-1,size(obj.pts_origin,1)) -...
                        circshift( eye(size(obj.pts_origin,1)-1,size(obj.pts_origin,1)), 1, 2);
                case 'reference'
                    obj.diff_mat = circshift( eye(size(obj.pts_origin,1)-1,size(obj.pts_origin,1)), 1, 2);
                    obj.diff_mat(:,1) = - 1;          
            end
            
            obj.n_channels = size(obj.diff_mat,1);
            
            obj.set_position([0,0,0], [0,0,0]);
            obj.set_linear_trajectory(0); %Defaul trajectory: no motion
        end

    end
    
end
    
