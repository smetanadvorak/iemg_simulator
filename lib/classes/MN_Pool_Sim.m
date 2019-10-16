classdef MN_Pool_Sim < handle
    % A class MotorNeuronPool_Sim is a container for the Motor Neuron Pool model
    
    properties
        motor_neurons;
        
        N ; %Number of mus
        rr; %Recruitment range: largest/smallest
        rm; %Recruitment maximum (when all the MUs are active)
        centers;
        
        sz; %Sizes
        rt; %Thresholds
        
        % Excitation-rate curve parameters
        type %Either linear rectified or deluca
        excfr_pars %For DeLuca model
        minfr; %Minimal firing rates
        maxfr; %Maximal firing rates
        frs;   %Firing rates slope  
        
        % IPI generation parameters
        CV = 1/6; % Coefficient of variation for gaussian model
        
        exc_fr_curves
        
    end
    
    methods
        function obj = MN_Pool_Sim(N, rr, rm)
            obj.N = N;
            obj.rr = rr;
            obj.rm = rm;
            obj.init_pool();
        end
        
        function init_pool(obj)
            [obj.sz, obj.rt] = generate_mu_recruitment_thresholds(obj.N, obj.rr, obj.rm);
            obj.sz = obj.sz(:);
            obj.rt = obj.rt(:);
        end
        
        function distribute_innervation_centers(obj, Rmuscle)
            n = 256;
            density_map = ones(n);
            [X,Y] = meshgrid(1:n,1:n);
            density_map(sqrt((X-n/2).^2 + (Y-n/2).^2) > n/2-1) = 0;
            
            vertices = zeros(2,obj.N);
            vertices(:,1) = [1;1];
            for i = 2:(obj.N+1)
                D = perform_fast_marching(1./density_map, vertices(:, 1:i-1));
                [~,ind] = max(D(:));
                [x,y] = ind2sub([n n],ind);
                vertices(:,i) = [x;y];
            end
            
            obj.centers = vertices(:,end:-1:2)';
            obj.centers = (obj.centers - n/2) / max(sqrt((obj.centers(:,1)-n/2).^2 + (obj.centers(:,2)-n/2).^2)) * Rmuscle;
        end
        
            
        function generate_minfr(obj, type, par1, par2)
            obj.type = type;
            obj.minfr = generate_mu_minfr( obj.rt/max(obj.rt), type, par1, par2 );
            obj.minfr = obj.minfr(:);
        end
        function generate_maxfr(obj, type, par1, par2)
            obj.type = type;
            obj.maxfr = generate_mu_maxfr( obj.rt/max(obj.rt), type, par1, par2 );
            obj.maxfr = obj.maxfr(:);
        end
        function generate_frs(obj, type, par1, par2)
            obj.type = type;
            obj.frs = generate_mu_fr_slope( obj.rt/max(obj.rt), type, par1, par2 );
            obj.frs = obj.frs(:);
        end
        
        function set_deluca_mdl(obj, type)
            % Alternative model of excitation-rate curves
            switch type
                case 'fdi'
                    A = 85;     B = 0.32;   C = -23;    D = 6.93;     E = 20.9;
                case 'vl'
                    A = 116;    B = 0.15;   C = -21;   D = 8.03;     E = 19.0;
                case 'ta'
                    A = 65;     B = 0.30;   C = -30;   D = 6.54;     E = 20.2;
            end
            obj.excfr_pars = [A; B; C; D; E];
            obj.type = 'deluca';
        end
        
        function fr = calculate_fr(obj, excitation, m)
            %fr = calculate_fr(obj.rt(k), obj.minfr(k), obj.maxfr(k), obj.frs(k), excitation);
            excitation = transpose(excitation(:));
            if nargin < 3
                if strcmp(obj.type, 'deluca')
                    dummy = num2cell(obj.excfr_pars);
                    [A,B,C,D,E] = deal(dummy{:});
                    fr = D*excitation + (C - A*exp(-excitation/B)).*obj.rt + E;
                    fr(excitation <= obj.rt) = 0;
                else
                    fr = obj.minfr + max(0, obj.frs .* (excitation - obj.rt));
                    fr = min(obj.maxfr, fr);
                    fr(excitation <= obj.rt) = 0;
                end
            else % Same but not vectorized
                if strcmp(obj.type, 'deluca')
                    dummy = num2cell(obj.excfr_pars);
                    [A,B,C,D,E] = deal(dummy{:});
                    fr = D*excitation + (C - A*exp(-excitation/B))*obj.rt(m) + E;
                    fr(excitation <= obj.rt(m)) = 0;
                else
                    fr = obj.minfr(m) + max(0, obj.frs(m) * (excitation - obj.rt(m)));
                    fr = min(obj.maxfr(m), fr);
                    fr(excitation <= obj.rt(m)) = 0;
                end
            end
            fr = transpose(fr); % Columns: excitation; Rows: motor neurons
        end

        %% Firings generation
        function [ spikes, next_state, last_ipi] = generate_spike_train_gauss(obj, T, prev_state, excitation, fs)
            spikes = zeros(numel(T),obj.N);
            ipi = zeros(obj.N,1);
            next_firing = prev_state;
            next_state = zeros(obj.N,1);
            last_ipi = zeros(obj.N,1);
            
            objrt = obj.rt; % Weird, but fast
            
            for m = 1:obj.N
                for t = 1:numel(T)
                    if excitation(t) > objrt(m)
                        if isnan(next_firing(m))
                            %ipi(m) = fs/calculate_fr(obj.rt(m), obj.minfr(m), obj.maxfr(m), obj.frs(m), excitation(t)); % mean IPI
                            ipi(m) = fs/obj.calculate_fr(excitation(t), m); % mean IPI
                            ipi(m) = ipi(m) + randn * ipi(m) * obj.CV; % add some variation
                            next_firing(m) = T(t) + round(ipi(m));
                        end
                        if T(t) == next_firing(m)
                            spikes(t,m) = 1;
                            %ipi(m) = fs/calculate_fr(obj.rt(m), obj.minfr(m), obj.maxfr(m), obj.frs(m), excitation(t));
                            ipi(m) = fs/obj.calculate_fr(excitation(t), m);
                            ipi(m) = ipi(m) + randn * ipi(m) * obj.CV;
                            next_firing(m) = T(t) + round(ipi(m));
                        end
                    else
                        next_firing(m) = nan;
                    end
                end
                last_ipi(m) = ipi(m);
                next_state(m) = next_firing(m);
            end
        end
        
        
        %% Vizualization methods
        function show_sizes(obj, ax)
            if nargin < 2
                ax = axes;
            end
            axes(ax); cla;
            
            stem(1:obj.N, obj.sz); hold on; 
            title('Thresholds and sizes distribution over the motor neuron pool');
            xlabel('Index of MU');
            ylabel('Motor neuron size, normalized units');
            yyaxis right;
            stem(1:obj.N,obj.rt, 'o');
            ylabel('Excitation rate, normalized units');
            legend({'MU sizes','Recruitment thresholds'});
            yyaxis left;
        end
        
        function show_centers(obj, ax, Rmuscle)
            if isempty(ax)
                ax = axes;
            end
            figure;
            axes(ax); cla;
            
            % Muscle border generation (should be passed as an argument in
            % future versions)
             phi_circle = linspace(0, 2*pi, 1000)';
             phi_circle = phi_circle(1:end-1); 
             muscle_border = [Rmuscle * cos(phi_circle), Rmuscle * sin(phi_circle)];
             plot(muscle_border(:,1), muscle_border(:,2), 'k', 'linewidth', 1); hold on;
            
            for i = 1:obj.N
                text(obj.centers(i,1)-0.1, obj.centers(i,2), num2str(i)); hold on;
                rad = sqrt(obj.sz(i)/sum(obj.sz) * pi * Rmuscle^2/pi);
                mu_area_circle = [rad * cos(phi_circle) + obj.centers(i,1), ...
                    rad * sin(phi_circle) + obj.centers(i,2)];
                
                plot(mu_area_circle(:,1), mu_area_circle(:,2), 'b', 'linewidth', 0.5);
            end
            axis equal
            axis tight
            %title('Motor neuron innervation centers and areas over the muscle cross-section');
            xlabel('x, mm'); ylabel('y, mm');
            %legend('Muscle border', 'Innervation areas');
        end
        
        function generate_exc_fr_curves(obj)
            excitation = linspace(0,1,1000)';
            obj.exc_fr_curves = zeros(length(excitation), obj.N);
            obj.exc_fr_curves = obj.calculate_fr(excitation);
            %for m = 1:obj.N
                %obj.exc_fr_curves(:,m) = obj.minfr(m) + max(0, obj.frs(m)*(exc - obj.rt(m)));
                %obj.exc_fr_curves(:,m) = min(obj.maxfr(m), obj.exc_fr_curves(:,m));
                %obj.exc_fr_curves(exc < obj.rt(m), m) = 0;
            %end
        end
        
        function pl = show_fr_exc_curves(obj, ax, inds)
            if nargin < 2 || isempty(ax)
                figure;
                ax = axes;
            end
            if nargin < 3 
                inds = 1:obj.N;
            end
            axes(ax); cla;
            
            obj.generate_exc_fr_curves();
            exc = linspace(0,1,size(obj.exc_fr_curves,1))';
            pl = plot(exc, obj.exc_fr_curves(:,inds), 'k', 'linewidth', 0.5);
            xlabel('Excitation, normalized'); ylabel('Firing rates, pulses per second');
        end
        
    end
    
end

