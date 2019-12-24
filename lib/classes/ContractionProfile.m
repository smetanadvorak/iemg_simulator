classdef ContractionProfile < handle
    %This class encapsulates the profile generation procedure, as well as
    %profile presentation
    
    properties
        fs      % Sampling frequency
        T       % Length in samples
        len     % Full length in seconds (including silent segements)
        act_len     % Length of the active part in seconds
        type
        
        profile  % Actual profile values
        
        timeline
        
        silent = nan;    % Silence time at start and end (in seconds)
        slope = nan;     % Slope value (%/s) for trapez and ramp profiles
        maxval = nan;    % Maximum value of profile
        minval = 0;      % Minimum value of profile (usually 0)
        
        recruitment_part; 
        derecruitment_part;
        constant_part;
        silent_part;
        
        fsl         % Subsampled frequency
        sub_profile % The same profile, subsampled
        sub_timeline;
    end
    
    methods
        function obj = ContractionProfile(len, fs, type, varargin)
            obj.fs = fs;
            obj.act_len = len;
            obj.type = type;
            
            switch obj.type
                case 'trapezoidal'
                    % Get the user's parameters of the profile
                    obj.maxval = varargin{1};
                    obj.minval = 0;
                    obj.slope = varargin{2}; % Slope of the trapezoidal curve, 1/s
                    
                    if numel(varargin) > 2
                        obj.silent = varargin{3};
                    else
                        obj.silent = 1; %1 second of silence in the beginning and in the end of profile
                    end
                    obj.len = obj.act_len + 2 * obj.silent;
                    
                    % Calculate the nodes of the profile
                    slope_len = obj.maxval/obj.slope;
                    plateau_len = max(0, obj.len - 2*(slope_len + obj.silent));
                    segments = cumsum([0, obj.silent, slope_len, plateau_len + eps, slope_len, obj.silent]);
                    
                    % Interpolate between nodes
                    obj.profile = [obj.minval, obj.minval, obj.maxval, obj.maxval, obj.minval, obj.minval];
                    [segments, inds] = unique(segments);
                    obj.profile = obj.profile(inds);
                    
                    obj.profile = interp1(segments * obj.fs, obj.profile, 1:round(obj.len*obj.fs));
                    obj.profile = obj.profile(:);
                    
                case 'constant'
                    % Get the user's parameters of the profile
                    obj.maxval = varargin{1};
                    obj.minval = 0;
                    if numel(varargin) > 1
                        obj.silent = varargin{2};
                    else
                        obj.silent = 1; %1 second of silence in the beginning 
                    end
                    obj.len = obj.act_len + obj.silent;
                    
                    segments = cumsum([0, obj.silent, obj.silent+eps, obj.len]);
                    obj.profile = [obj.minval, obj.minval, obj.maxval, obj.maxval];
                    obj.profile = interp1(segments * obj.fs, obj.profile, 1:round(obj.len*obj.fs));
                    obj.profile = obj.profile(:);
                    
                case 'ramp'
                    % Get the user's parameters of the profile
                    obj.maxval = varargin{1};
                    obj.minval = 0;
                    obj.slope = varargin{2}; % Slope of the curve, %/s
                    
                    if numel(varargin) > 2
                        obj.silent = varargin{3};
                    else
                        obj.silent = 1; %1 second of silence in the beginning
                    end
                    obj.len = obj.act_len + obj.silent;
                    
                    if strcmp(obj.slope, 'auto')
                        obj.slope = (obj.maxval - obj.minval) / (obj.len - obj.silent);
                    end
                    
                    % Calculate the nodes of the profile
                    slope_len = obj.maxval/obj.slope;
                    if slope_len + obj.silent > obj.len
                        warning('Ramp profile slope, max and length are specified so that the maximum value is not reached. Consider setting longer profile duration or smaller slope.');
                        slope_len = len - obj.silent;
                    end
                    
                    plateau_len = max(0, obj.len - (slope_len + obj.silent));
                    segments = cumsum([0, obj.silent, slope_len, plateau_len+eps]);
                    
                    % Interpolate between nodes
                    obj.profile = [obj.minval, obj.minval, obj.maxval, obj.maxval];
                    [segments, inds] = unique(segments);
                    obj.profile = obj.profile(inds);
                    obj.profile = interp1(segments * obj.fs, obj.profile, 1:round(obj.len*obj.fs));
                    obj.profile = obj.profile(:);
                    
                case "rect"
                    % Get the user's parameters of the profile
                    obj.maxval = varargin{1};
                    obj.minval = 0;
                    
                    if numel(varargin) > 1
                        obj.silent = varargin{2};
                    else
                        obj.silent = 1; %1 second of silence in the beginning and in the end of profile
                    end
                    obj.len = obj.act_len + 2 * obj.silent;
                    
                    % Calculate the nodes of the profile
                    plateau_len = max(0, obj.len - 2*obj.silent);
                    segments = cumsum([0, obj.silent, eps, plateau_len, eps, obj.silent]);
                    
                    % Interpolate between nodes
                    obj.profile = [obj.minval, obj.minval, obj.maxval, obj.maxval, obj.minval, obj.minval];
                    [segments, inds] = unique(segments);
                    obj.profile = obj.profile(inds);
                    
                    obj.profile = interp1(segments * obj.fs, obj.profile, 1:round(obj.len*obj.fs));
                    obj.profile = obj.profile(:);
                    
                    
                %case 'free'
                
%                 case 'sinusoidal'
%                     % Get the user's parameters of the profile
%                     meanval = varargin{1};
%                     ampl = varargin{2};
%                     
%                     obj.maxval = meanval+ampl;
%                     obj.minval = 0;
%                     
%                     if numel(varargin) > 2
%                         obj.silent = varargin{2};
%                     else
%                         obj.silent = 1; %1 second of silence in the beginning and in the end of profile
%                     end
%                     obj.len = obj.act_len + 2 * obj.silent;
                   
                    
                otherwise
                    error('Profile generation failed because specified profile isn''t recognized');
            end
               
            obj.timeline = (1:length(obj.profile))/obj.fs;
            obj.T = length(obj.timeline); %Length in samples
            
            obj.recruitment_part = diff([obj.profile(1,:); obj.profile]) > 0;
            obj.derecruitment_part = diff([obj.profile(1,:); obj.profile]) < 0;
            obj.silent_part = obj.profile == 0;
            obj.constant_part = diff([obj.profile(1,:); obj.profile]) == 0;
            obj.constant_part(obj.silent_part) = 0;
        end
        
        function subsample(obj, fsl)
            obj.fsl = fsl;
            obj.sub_profile = obj.profile((1 : round(obj.T/obj.fs*obj.fsl)) * round(obj.fs/obj.fsl));
            obj.sub_timeline = (1 : round(obj.T/obj.fs*obj.fsl)) * round(obj.fs/obj.fsl) / obj.fs;
        end
        
        function str = summary(obj)
            str = ['Type: ', obj.type, ', length: ', num2str(obj.len-obj.silent), ' sec, maximum value: ', num2str(obj.maxval)];
        end
        
        function show(obj, ax)
            if nargin < 2 || isempty(ax)
                figure;
                ax = axes;
            end
            %plot(ax, obj.timeline, obj.profile, 'b', 'linewidth',2);
            
            plot(ax, obj.timeline(obj.recruitment_part), obj.profile(obj.recruitment_part), 'g.', 'markersize',5, 'DisplayName', 'Recruitment segment'); hold on;
            plot(ax, obj.timeline(obj.derecruitment_part), obj.profile(obj.derecruitment_part), 'r.', 'markersize',5, 'DisplayName', 'Derecruitment segment');
            plot(ax, obj.timeline(obj.silent_part), obj.profile(obj.silent_part), 'k.', 'markersize',5, 'DisplayName', 'Silent segment');
            plot(ax, obj.timeline(obj.constant_part), obj.profile(obj.constant_part), 'b.', 'markersize',5, 'DisplayName', 'Constant segment');
            axis tight
            title('Contraction profile');
            ylabel('Level of contraction, per cent MVC');
            xlabel('Time, s');
            legend(gca, 'show');%'Recruitment segment', 'Derecrtuitment segement', 'Silent segment', 'Constant segment');            
        end
        
    end
end

