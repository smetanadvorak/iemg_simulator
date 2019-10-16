% ToDo: all the estimators should downsample the data themselves. Or you do
% the feature extracture (activity, firing rates) outside of them (like
% now)

%% Define the model type, configure it
fe_training_lag = round(75/1000*fsl);
fe_inference_lag = round(1/1000*fsl);
fe_window = ones(ceil(250/1000*fsl), 1); 
%fe_window = gausswin(ceil(250/1000*fsl), 1); 
fe_window = fe_window./sum(fe_window);

%% Define effect estimator
%estimator = PC_CDR_online(fe_window, fe_training_lag, fe_inference_lag);
%estimator = PC_RS_online(fe_training_lag, fe_inference_lag); 
estimator.dist_type = 'logistic';
%estimator = PC_RC_online('deluca', fe_training_lag, fe_inference_lag);
estimator = PC_RSRC_online(fe_training_lag, fe_inference_lag);
clear fr_estimator


%% Define firing rate estimators
for m = 1:size(train_spk,2)
    fr_estimator(m) = IFR_online(fe_window);
    %fr_estimator(m) = IFR_online_weibull(fe_window, 30/1000 * fsl, 100/1000 * fsl, 10, 1, 1);
end

%% Train the estimator
set_exact_estimator_parameters = 1;
if set_exact_estimator_parameters
    estimator.set_exact_parameters(mu_pool.mn_pool, mf_mdl, detectable_ind);
    warning('Estimator''s parameters are manually set to exact values!');
    
else % Do training on training dataset    
    % Calculate firing rates
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
    
    
    %% Predict training data
    train_predicted_dwns = zeros(size(train_aux_dwns));
    switch estimator.type
        case 'RS'
            for i = 1:numel(train_predicted_dwns)
                estimator.update(train_act_dwns(i,:));
                train_predicted_dwns(i) = estimator.predict();
            end
            
        case 'RC'
            for i = 1:numel(train_predicted_dwns)
                estimator.update(train_fr_dwns(i,:));
                train_predicted_dwns(i) = estimator.predict();
            end
            
        case 'RSRC'
            for i = 1:numel(train_predicted_dwns)
                estimator.update(train_fr_dwns(i,:), train_act_dwns(i,:));
                train_predicted_dwns(i) = estimator.predict();
            end
            
        case 'CDR'
            for i = 1:numel(train_predicted_dwns)
                estimator.update(train_spk_dwns(i,:));
                train_predicted_dwns(i) = estimator.predict();
            end
    end

    train_predicted = interp1(train_predicted_dwns, linspace(0, length(train_predicted_dwns), length(train_aux))');
    fe_R2_train = 1-var(train_aux - train_predicted)/var(train_aux);
    %% Plot training data
    figure;
    train_timeline = (1:length(train_aux_dwns))/fs;
    plot(train_timeline, train_aux_dwns, 'k', 'linewidth',2); hold on;
    plot(train_timeline, train_predicted_dwns, 'r');
    ylabel('Force, normalized');
    yyaxis right;
    ylim([0,inf]);
    %plot(train_timeline(:)', train_fr_dwns, '-');
    xlabel('Time, s'); ylabel('Firing rate, Hz'); axis tight
end
