cd simulation_output
cd(output_dir_name);

fid = fopen('simulation_report.txt', 'w');


% Preambule
fprintf(fid, 'This EMG signal is simuluated using simulator created by Konstantin Akhmadeev.\n\n');


% Motor neuron pool:
fprintf(fid, 'Number of Motor Neurons: %d\n', N);
fprintf(fid, 'Recruitment range: %d\n', rr);
fprintf(fid, 'Motor neuron pool Model: see Fuglevand, Winter, Patla - 1993\n\n');


% Muscle parameters:
fprintf(fid, 'Muscle form: cylindrical\n');
fprintf(fid, 'Muscle radius: %d mm\n', Rmuscle);
fprintf(fid, 'Muscle length: %d mm\n', Lmuscle);
fprintf(fid, 'Number of muscle fibers: %d\n', Nmf);

mf_diameters = []; mf_cv = [];
for i = 1:N, mf_diameters = [mf_diameters; MUs(i).mf_diameters]; end;
for i = 1:N, mf_cv = [mf_cv; MUs(i).mf_cv]; end;
fprintf(fid, 'Muscle diameters distribution: mu = %2.2f mu, std = %2.2f mu\n', mean(mf_diameters)*1000, std(mf_diameters)*1000);
fprintf(fid, 'Conduction velocity distribution: mu = %2.2f m/s, std = %2.2f m/s\n\n', mean(mf_cv)/1000, std(mf_cv)/1000);


% Neuromuscular junction parameters:
fprintf(fid, 'Neuromuscular jitter std: %2.2f\n\n', nmj_jitter);


% Electrode parameters: 
fprintf(fid, 'Electrode type: %s\n', electrode_type);
fprintf(fid, 'Number of points: %d\n', size(electrode_pts_init, 1));
fprintf(fid, 'Number of trajectory nodes: %d\n\n', traj_n_nodes);


% Profile parameters:
fprintf(fid, 'Profile type: %s\n', profile_type);
fprintf(fid, 'Percent of MVC: %d\n', percent_mvc);
fprintf(fid, 'Profile duration: %d sec\n\n', round(T/fs));


% PID parameters:


% EMG signal parameters:
fprintf(fid, 'EMG SNR (relative to MVC EMG): %d dB\n\n', SNR);

% Annotation parameters:
fprintf(fid, 'Number of active detectable MUs: %d \n\n', sum(any(spikes(:,detectable_ind))));


cd ../..