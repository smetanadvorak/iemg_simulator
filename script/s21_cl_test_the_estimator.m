%% If needed, set exact simulation parameters to the estimator
plot_criterion = 1;

%% Downsample the testing data
downsampler_aux.reset();
downsampler_spk.reset();
downsampler_act.reset();

test_aux_dwns = zeros(size(test_aux));
test_act_dwns = zeros(size(test_spk));
test_spk_dwns  = zeros(size(test_spk));
test_act = spikes2activity(test_spk, length(fe_window)/fsl*fs, 1);

for t = 1:size(test_spk,1)
    test_aux_dwns(t) = downsampler_aux.update(test_aux(t));
    test_act_dwns(t,:) = downsampler_act.update(test_act(t,:));
    test_spk_dwns(t,:) = downsampler_spk.update(test_spk(t,:));
end

test_aux_dwns(isnan(test_aux_dwns(:,1)),:) = [];
test_act_dwns(isnan(test_act_dwns(:,1)),:) = [];
test_spk_dwns(isnan(test_spk_dwns(:,1)),:) = [];

test_act_dwns = double(test_act_dwns > 0);
test_spk_dwns = double(test_spk_dwns > 0);

% Calculate firing rates
test_fr_dwns = zeros(size(test_spk_dwns));
for m = 1:size(test_spk_dwns,2)
    for t = 1:size(test_spk_dwns,1)
        test_fr_dwns(t,m) = fr_estimator(m).update(test_spk_dwns(t,m)) * fsl;
    end
    fr_estimator(m).reset();
end

%% Allocate crit matrices
test_predicted_dwns = zeros(size(test_aux_dwns));
crit_nodes = linspace(0, max(test_aux_dwns)*1.25 ,500);
crit_hist = zeros(length(test_aux_dwns), numel(crit_nodes));
crit_low = zeros(length(test_aux_dwns), numel(crit_nodes));
crit_top = zeros(length(test_aux_dwns), numel(crit_nodes));
crit_gauss = zeros(length(test_aux_dwns), numel(crit_nodes));

%% Test the estimator 

estimator.reset();
fprintf('\nTesting the estimator...');

switch estimator.type
    case 'RS'
        for i = 1:numel(test_predicted_dwns)
            estimator.update(test_act_dwns(i,:));
            test_predicted_dwns(i) = estimator.predict();
            [crit_hist(i, :), crit_low(i,:), crit_top(i,:)] = estimator.get_inference_criterion(test_act_dwns(i,:), crit_nodes);
        end
        
    case 'RC'
        for i = 1:numel(test_predicted_dwns)
            estimator.update(test_fr_dwns(i,:));
            test_predicted_dwns(i) = estimator.predict();
            [crit_hist(i, :)] = estimator.get_inference_criterion(test_fr_dwns(i,:), test_act_dwns(i,:), crit_nodes);
        end
        
    case 'RSRC'
        for i = 1:numel(test_predicted_dwns)
            estimator.update(test_fr_dwns(i,:), test_act_dwns(i,:));
            test_predicted_dwns(i) = estimator.predict();
            [crit_hist(i, :), crit_low(i,:), crit_top(i,:), crit_gauss(i,:)] = estimator.get_inference_criterion(test_fr_dwns(i,:), test_act_dwns(i,:), crit_nodes);
        end
        
    case 'CDR'
        for i = 1:numel(test_predicted_dwns)
            estimator.update(test_spk_dwns(i,:));
            test_predicted_dwns(i) = estimator.predict();
        end        
end

test_predicted = interp1(test_predicted_dwns, linspace(0, length(test_predicted_dwns), length(test_aux))');
fe_R2_test = 1-var(test_aux - test_predicted)/var(test_aux);
fprintf(' Done!\n');

%% Plot the test results
figure; 
%plot(test_timeline_dwns, test_fr_dwns, '-');
%legend('Observation', 'Predicted', 'Rates');



% Plot firing rates
test_timeline = (1:length(test_aux))/fs;
test_timeline_dwns = (1:length(test_aux_dwns))/fsl;
plot(test_timeline, test_aux, 'k', 'linewidth',2); hold on; 
plot(test_timeline, test_predicted,'g', 'linewidth',2);
xlabel('Time, s'); ylabel('Force, normalized'); axis tight
legend('True', 'Estimated');

% Plot spikes
yyaxis right
plot(test_timeline, test_spk .* repmat(1:size(test_spk,2), size(test_spk,1), 1), '.')
%plot(test_timeline, test_spk .* repmat(estimator.pars(:,1)', size(test_spk,1), 1), '.')
ylim([0.5, 0.5 + find(any(test_spk), 1, 'last')]);
ylabel('Motor neuron');
%set(gca, 'ytick', []);
legend('Observation', 'Predicted', 'Spikes');


%% Plot criterion
%crit_inds = conv(sum(test_act_dwns, 2) == 0, ones(51,1), 'same') > 0;
if plot_criterion
    figure; 
    crit_hist_viz = -log(1+crit_hist);
    crit_hist_viz(any(diff(test_act_dwns),2),:)=nan;
    % Plot criterion
    [X,Y] = meshgrid((1:length(test_aux_dwns))'/fsl, crit_nodes);
    mesh(X, Y, crit_hist_viz', 'facealpha',0.85); hold on
    
    % Plot criterion maximal points
    [crit_mins, min_inds] = max(crit_hist_viz,[],2);
    %plot3((1:length(test_aux_dwns))'/fsl, crit_nodes(min_inds), crit_mins, 'r.', 'markersize', 15);
    plot3(test_timeline, test_predicted, 0.15*ones(size(test_timeline)), 'linewidth',3, 'color', 'g');
    plot3(test_timeline, test_aux, 0.05*ones(size(test_timeline)), 'linewidth',3, 'color', 'k');
    xlabel('Time, s'); ylabel('Effect $e$, normalized'); zlabel('$-\log C_I$'); 
end

clear X Y

%%
%clear signals_file annotation_file fe_* clear EMGSIGNAL AUXSIGNAL current_path
%clear test_timeline test_aux* test_predicted* test_activity* i
%clear train_* test_*