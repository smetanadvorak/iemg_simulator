
%% Define the model type, configure it
fe_training_lag = round(100/1000*fs);
fe_inference_lag = round(0/1000*fs);
fe_smoothing_window = gausswin(200/1000*fs);
fe_smoothing_window = fe_smoothing_window./sum(fe_smoothing_window);
estimator = PC_CDR_online(fe_smoothing_window, fe_training_lag, fe_inference_lag); 


%% Train
fprintf('\nTraining the estimator...');
estimator.train(train_spikes, train_aux);
fprintf(' Done!\n');


%% Test the estimator 
fprintf('\nTesting the estimator...');
test_predicted = zeros(size(test_aux));
estimator.reset();
for i = 1:numel(test_predicted)
    estimator.update(test_spikes(i,:));
    test_predicted(i) = estimator.predict();
end

fe_R2_test = 1-var(test_aux - test_predicted)/var(test_aux);
fprintf(' Done!\n');



%% Plot the test results
figure; 
test_timeline = (1:length(test_aux))/fs;
plot(test_timeline, test_aux, 'k', 'linewidth',2); hold on; 
plot(test_timeline, test_predicted,'g');
xlabel('Time, s'); ylabel('Force, normalized'); axis tight
legend('Observation', 'Predicted (CDR)');
title(sprintf('Force estimation from CDR, $R^2$=%2.3f. Test data.', fe_R2_test)); figure2page('A4')


%%
clear signals_file annotation_file fe_* clear EMGSIGNAL AUXSIGNAL current_path