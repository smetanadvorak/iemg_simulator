fsd = 25000; % Downsampling frequency

ROI = [1; round(T/fs)-1]; % not implemented


if fsd < fs
    %% Resample EMG
    emg_downsample = EMG_filtering_fir(emg, [], fsd/2, [], fs);
    emg_downsample = resample(emg_downsample, fsd, fs);
    
    emg_clear_downsample = EMG_filtering_fir(emg_clear, [], fsd/2, [], fs);
    emg_clear_downsample = resample(emg_clear_downsample, fsd, fs);
    
    %% Resample Force
    aux_downsample = EMG_filtering_fir(force, [], fsd/4, [], fs);
    aux_downsample = resample(aux_downsample, fsd, fs);
    
    %% Resample Excitation
    excitation_downsample = EMG_filtering_fir(excitation, [], fsd/4, [], fs);
    excitation_downsample = resample(excitation_downsample, fsd, fs);
    
    %% Resample the profile
    profile = resample(profile, fsd, fs);
    
    %% Resample Annotation
    % No need: ann format keeps true times.
    
    %% Resample Dictionary
%     dictionary_full_init_downsample = zeros();
%     for m = 1:N
%         dictionary_full_init_downsample(:,:,m) = resample(dictionary_full(:,:,m), fsd, fs);
%     end
%     dictionary_detectable_downsample = zeros(size(dictionary_detectable));
%     for m = 1:size(dictionary_detectable_downsample,2)
%         dictionary_detectable_downsample(:,:,m) = resample(dictionary_detectable(:,:,m), fsd, fs);
%     end
    
    
else
    excitation_downsample = excitation;
    emg_downsample = emg;
    emg_clear_downsample = emg_clear;
    aux_downsample = force;
    profile_downsample = profile;
end

%% Save stuff
%cd ..

%mkdir 'simulation_output';
cd 'simulation_output';

output_dir_name = [profile_type, '_', electrode_type_short, '_', num2str(N), '_MUs'];

mkdir(output_dir_name);
cd(output_dir_name);


% ================= Save EMG and AUX signals =================
EMGSIGNAL.data = emg_downsample; %Well... fake amplitude
EMGSIGNAL.rate = fsd;

AUXSIGNAL.data = aux_downsample;
AUXSIGNAL.rate = fsd;

EMGSIGNAL_clear.data = emg_clear_downsample;
EMGSIGNAL_clear.rate = fsd;

EXCSIGNAL.data = excitation_downsample;
EXCSIGNAL.rate = fsd;
PRFSIGNAL.units = 'Normalized to maximum value';

PRFSIGNAL.data = profile_downsample;
PRFSIGNAL.rate = fsd;
PRFSIGNAL.units = 'Normalized to MVC';
PRFSIGNAL.type = profile_type;



save('simulation_signals', 'EMGSIGNAL', 'AUXSIGNAL', 'EMGSIGNAL_clear', ...
     'EXCSIGNAL', 'PRFSIGNAL', 'mvc_emg_std', 'detectable_ind');
%save('simulation_muscle', 'sz', 'mf*', 'electrode*', 'MUs', 'force_profile');
%    save('sim_dictionary_full', 'dictionary_full_downsample');
%    save('sim_dictionary_detectable', 'dictionary_detectable_downsample');


% ================= Save dictionnaries ========================
save('dictionary_full', 'dictionary_full_init', 'comment_full_init', ...
                        'dictionary_full_traj', 'comment_full_traj');
                    
save('dictionary_detectable', 'dictionary_detectable_init', 'comment_detectable_init',...
                              'dictionary_detectable_traj', 'comment_detectable_traj');

% ================= Save annotation (nmj time) ===========================
%firings2ann(annotation_causal_full, 'sim_annotation_causal_full.ann', fs);
%firings2ann(annotation_causal_detectable, 'sim_annotation_causal_detectable.ann', fs);

% ================= Save annotation (ap center time) ============================
firings2ann(firings_centered_full, 'sim_annotation_centered_full', fs);
firings2ann(firings_centered_detectable, 'sim_annotation_centered_detectable', fs);
    
cd ../..
%cd 'simulator'
clear emg_downsample aux_downsample