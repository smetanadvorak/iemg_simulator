% ToDo: all the estimators should downsample the data themselves. Or you do
% the feature extracture (activity, firing rates) outside of them (like
% now)

%% Use training data or set exact parameters to the estimator
set_exact_estimator_parameters = 0;

%% Define the model type, configure it
fe_training_lag = round(75/1000*fsl);
fe_inference_lag = round(1/1000*fsl);
fe_window = ones(250/1000*fsl, 1); 
fe_window = fe_window./sum(fe_window);

%% Define effect estimator
%estimator = PC_CDR_online(fe_window, fe_training_lag, fe_inference_lag);
estimator = PC_RS_online(fe_training_lag, fe_inference_lag);
%estimator = PC_RC_online('linear', fe_training_lag, fe_inference_lag);
%estimator = PC_RSRC_online(fe_training_lag, fe_inference_lag);
clear fr_estimator


%% Define firing rate estimators
for m = 1:size(train_spikes,2)
    fr_estimator(m) = IFR_online(fe_window);
    %fr_estimator(m) = IFR_online_weibull(fe_window, 30/1000 * fsl, 100/1000 * fsl, 10, 1, 1);
end

%% Train the estimator
if set_exact_estimator_parameters
    estimator.set_exact_parameters(mu_pool.mn_pool, mf_mdl, detectable_ind);
    warning('Estimator''s parameters are manually set to exact values!');
    
else % Do training on training dataset    
    %% Calculate firing rates
    train_fr_dwns = zeros(size(train_spk_dwns));
    for m = 1:size(train_spk_dwns,2)
        for t = 1:size(train_spk_dwns,1)
            train_fr_dwns(t,m) = fr_estimator(m).update(train_spk_dwns(t,m)) * fsl;
        end
        fr_estimator(m).reset();
    end
    
    %% Train
    fprintf('\nTraining the estimator...');
    switch estimator.type
        case 'RS'
            estimator.train(train_act_dwns, train_aux_dwns);
        
        case 'RC'
            estimator.train(train_fr_dwns, train_aux_dwns);
            
        case 'RSRC'
            estimator.train(train_fr_dwns, train_act_dwns, train_aux_dwns);
            
        case 'CDR'
            estimator.train(train_spk_dwns, train_aux_dwns);
    end
    fprintf(' Done!\n');
    
    %% Plot training data
    figure;
    train_timeline = (1:length(train_aux_dwns))/fs;
    plot(train_timeline, train_aux_dwns, 'k', 'linewidth',2); hold on;
    ylabel('Force, normalized');
    yyaxis right;
    plot(train_timeline(:)', train_fr_dwns, '-');
    xlabel('Time, s'); ylabel('Firing rate, Hz'); axis tight
end
