classdef MuscleForceMdl < handle
    % This class encapsulates the:
    % - twitch models for N MUs
    % - faster-than-linear coefficients for force calculation
    % - offline force calculation
    % - online force calculation
   
    
    properties
        % Twitches parameters
        Tmax;           % Maximum time-to-peak delay (see Fuglevand)
        Tr = 3;         % Time to peak range
        Tcoeff;
        T;              % Times to peak
        P;              % Peak heights
        
        % Twitches
        twitch_mat;
        twitch_cell;
        
        % General
        hfs;            % Modelling frequency (high)
        N;              % Number of motor neurons
        size_magnitude
        fmax;           % Maximal voluntary contraction (for normalization)
        online_buffer;
        
        % Quasistatic exc2force model
        qs_e2f_mdl;
        % Quasistatic force2exc model
        qs_f2e_mdl;
    end
    
    methods
        function obj = MuscleForceMdl(N, size_magnitude, hfs)
            obj.hfs = hfs;
            obj.N = N;
            obj.size_magnitude = size_magnitude;
            obj.init_twitch_parameters();
            obj.init_twitches();
        end
        
        function init_twitch_parameters(obj)
            %(see Fuglevand et al. - Models of Recruitment and Rate coding organisation in motor units, 1993.)
            % Peak twitch forces: proportional to sizes
            % Fuglevand's model:
            %           twitch(t) = Pf/Tf * t .* exp(1 - t/Tf);
            % Tf is time to peak;
            % Pf is peak height;
            
            % Peak heights:
            obj.P = exp((log(obj.size_magnitude)/obj.N) * (1:obj.N));
            
            % Times to peak:
            obj.Tmax = 90/1000 * obj.hfs; % Maximum time-to-peak delay: 90 ms (Fuglevand)
            obj.Tr = 3; % Time to peak range
            obj.Tcoeff = log(obj.size_magnitude)/log(obj.Tr);
            obj.T = obj.Tmax * (1./obj.P).^(1/obj.Tcoeff); % see Eq. 15 of Fuglevand
        end
        
        
        function init_twitches(obj)
            % Precaclulate the twitches and save them to a matrix(equalized
            % lengths) and cell (individual effective lengths).
            max_twitch_len = ceil(5*max(obj.T));
            twitch_timeline = (0:(max_twitch_len-1))';
            obj.twitch_cell = cell(obj.N,1);
            obj.twitch_mat = zeros(max_twitch_len, obj.N);
            
            for i = 1:obj.N
                twitch = twitch_fuglevand(twitch_timeline, obj.P(i), obj.T(i));
                obj.twitch_mat(:,i) = twitch;
                
                effective_length = min(length(twitch), ceil(5*obj.T(i)));
                obj.twitch_cell{i} = twitch(1:effective_length);
            end
        end
        
        
        function twitch = get_twitch(~, t, P, T)
            twitch = P/T * t .* exp(1 - t/T); 
        end
        
        function gain = get_gain(~, ipi, T)
            % gain = get_gain(ipi, T) returns the gain value for the
            % force output for a motor unit with current inter-pulse-interval ipi and
            % T-parameter of the twitch T. This function corresponds to Fuglevands
            % nonlinear gain model for the force output, see Fuglevant - Models of Rate
            % Coding ..., eq. 17.
            
            Sf = @(x)(1 - exp( -2*(x).^3));
            
            inst_dr = T./ipi;           % Instantaneous discharge rate
            gain = ones(size(inst_dr)); % Gain
            gain(inst_dr > 0.4) = Sf(inst_dr(inst_dr > 0.4)) ...
                                   ./ inst_dr(inst_dr > 0.4) ...
                                    / (Sf(0.4)/0.4);  % See Fuglevand, eq. 17
        end
        
        
        function hf = plot_twitches(obj, ax)
            if nargin < 2
                figure;
                ax = axes;
            end
            timeline = (1:size(obj.twitch_mat,1))/obj.hfs;
            plot(ax, timeline, obj.twitch_mat, 'k');
            title("Twitch waveforms"); xlabel("Time, s"); ylabel("Amplitude");
            hf = gcf;
        end
        
        %%
        function force_hfs = generate_force_offline(obj, spikes, prefix)
            if nargin < 3
                prefix = '';
            end
            
            L = size(spikes,1);
            
            % IPI signal generation out of spikes signal (for gain nonlinearity)
            [~, ipi] = sawtooth2ipi(spikes2sawtooth([spikes(2:end, :); zeros(1,obj.N)]));
            gain = nan(size(spikes));
            for n = 1:obj.N
                gain(:,n) = obj.get_gain(ipi(:,n), obj.T(n));
            end
            
            % Generate force
            force_hfs = zeros(L,1);
            for n = 1:obj.N
                for t = 1:L
                    if spikes(t,n)
                        twitch_to_add = obj.twitch_cell{n};
                        to_take = min(length(twitch_to_add), L-t+1);
                        force_hfs(t:(t+to_take-1))  = ...
                            force_hfs(t:(t+to_take-1)) + gain(t,n) * twitch_to_add(1:to_take);
                    end
                end
                fprintf('%s %d Twitch trains are generated\n',prefix, n);
            end
            force_hfs = force_hfs/obj.fmax; % Normalization to % MVC
        end
        
        
        %%
        function init_online_buffer(obj, len)
            if nargin < 2
                len = size(obj.twitch_mat, 1);
            end
            obj.online_buffer = zeros(len, 1);
        end
        
        %%
        function force = generate_force_online(obj, spikes, ipi)
            if ~any(size(spikes) == 1) || ~any(size(ipi) == 1)
                error('To work in online manner, pass the spikes vector element by element');
            end
            spikes = spikes(:); ipi = ipi(:);
            
            for i = 1:obj.N
                if spikes(i)
                    to_add = obj.get_gain(ipi(i), obj.T(i))   *   obj.twitch_mat(:,i);
                    obj.online_buffer = obj.online_buffer + to_add;
                end
            end
            
            % Get the current contraction
            force = obj.online_buffer(1) / obj.fmax;
            % Shift buffer
            %obj.online_buffer(1) = []; obj.online_buffer(end+1) = 0;
            obj.online_buffer = [obj.online_buffer(2:end); 0];
        end
        
        %%
        function [mvc_force, fmax] = normalize_mvc(obj, spikes)
            obj.fmax = 1;
            try
            mvc_force = obj.generate_force_offline(spikes, 'MVC measurement:');
            fmax = mean(mvc_force(round(end/2):end));
            catch err
                obj.fmax = [];
                rethrow(err);
            end
            obj.fmax = fmax;
            mvc_force = mvc_force/obj.fmax;
        end
        
        
        function init_quasistatic_e2f_f2e_models(obj, mu_pool)
            qsi_T = 25*obj.hfs;
            qsi_excitation = linspace(0,1,qsi_T)';
            qsi_spikes = mu_pool.mn_pool.generate_spike_train_gauss(1:qsi_T, nan(mu_pool.mn_pool.N,1), qsi_excitation, obj.hfs);
            qsi_force = obj.generate_force_offline(qsi_spikes);
            
            % Invert the curve, get weightened polynomial interpolation that passes through zero
            w = 1.25-qsi_excitation;

            obj.qs_f2e_mdl = (w.*[qsi_force.^5, qsi_force.^4, qsi_force.^3, qsi_force.^2, qsi_force.^1])\(w.*qsi_excitation);
            obj.qs_f2e_mdl = [obj.qs_f2e_mdl; 0];
            
            obj.qs_e2f_mdl = (w.*[qsi_excitation.^5, qsi_excitation.^4, qsi_excitation.^3, qsi_excitation.^2, qsi_excitation.^1])\(w.*qsi_force);
            obj.qs_e2f_mdl = [obj.qs_e2f_mdl; 0];
        end
        
        function res = e2f(obj, e)
            res = polyval(obj.qs_e2f_mdl, e);
        end
        
        function res = f2e(obj, f)
            res = polyval(obj.qs_f2e_mdl, f);
        end
            
    end
end

