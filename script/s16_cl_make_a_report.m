cd simulation_output
cd(output_dir_name);

fid = fopen('simulation_report.txt', 'w');


% Preambule
fprintf(fid, 'This EMG signal is simuluated using simulator created by Konstantin Akhmadeev.\n\n');


% Motor neuron pool:
fprintf(fid, 'Number of Motor Neurons: %d\n', mu_pool.mn_pool.N);
fprintf(fid, 'Recruitment range: %d\n', mu_pool.mn_pool.rr);
fprintf(fid, 'Motor neuron pool Model: see Fuglevand, Winter, Patla - 1993\n\n');

% Muscle parameters:
fprintf(fid, 'Muscle form: cylindrical\n');
fprintf(fid, 'Muscle radius: %d mm\n', Rmuscle);
fprintf(fid, 'Muscle length: %d mm\n', Lmuscle);
fprintf(fid, 'Number of muscle fibers: %d\n', mu_pool.Nmf);

fprintf(fid, 'Muscle diameters distribution: mu = %2.2f mu, std = %2.2f mu\n', mean(mu_pool.mf_diameters)*1000, std(mu_pool.mf_diameters)*1000);
fprintf(fid, 'Conduction velocity distribution: mu = %2.2f m/s, std = %2.2f m/s\n\n', mean(mu_pool.mf_cv)/1000, std(mu_pool.mf_cv)/1000);


% Neuromuscular junction parameters:
fprintf(fid, 'Neuromuscular jitter std: %2.2f\n\n', nmj_jitter);


% Electrode parameters: 
fprintf(fid, 'Electrode type: %s\n', electrode.type);
fprintf(fid, 'Number of points: %d\n', size(electrode.pts_init, 1));
fprintf(fid, 'Number of trajectory nodes: %d\n\n', electrode.n_nodes);


% Profile parameters:
fprintf(fid, 'Profile type: %s\n', profile.type);
fprintf(fid, 'Percent of MVC: %d\n', profile.maxval * 100);
fprintf(fid, 'Profile duration: %d sec\n\n', round(T/fs));


% PID parameters:


% EMG signal parameters:
fprintf(fid, 'EMG SNR (relative to MVC EMG): %d dB\n\n', SNR);

% Annotation parameters:
fprintf(fid, 'Number of active detectable MUs: %d \n\n', sum(any(spikes(:,detectable_ind))));


cd ../..