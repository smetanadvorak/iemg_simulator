%% If needed, set exact simulation parameters to the estimator
plot_criterion = 1;

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

test_predicted = interp1(test_predicted_dwns, linspace(0,length(test_predicted_dwns), length(test_predicted_dwns)/fsl*fs)');
fe_R2_test = 1-var(test_aux - test_predicted)/var(test_aux);
fprintf(' Done!\n');

%% Plot the test results
figure; 
test_timeline = (1:length(test_aux))/fs;
test_timeline_dwns = (1:length(test_aux_dwns))/fsl;
plot(test_timeline, test_aux, 'k', 'linewidth',2); hold on; 
plot(test_timeline, test_predicted,'g', 'linewidth',2);
xlabel('Time, s'); ylabel('Force, normalized'); axis tight

% Plot spikes
% yyaxis right
% plot(test_timeline, test_spikes .* repmat(1:size(test_spikes,2), size(test_spikes,1), 1), 'o')
% ylim([0.5, 0.5 + find(any(test_spikes), 1, 'last')]);
% set(gca, 'ytick', []);
% legend('Observation', 'Predicted', 'Spikes');

% Plot firing rates
yyaxis right
plot(test_timeline_dwns, test_fr_dwns, '-');
legend('Observation', 'Predicted', 'Rates');

%% Plot criterion
%crit_inds = conv(sum(test_act_dwns, 2) == 0, ones(51,1), 'same') > 0;
if plot_criterion
    figure; 
    crit_hist_viz = -log(crit_hist);
    
    % Plot criterion
    [X,Y] = meshgrid(crit_nodes(crit_nodes_inds), (1:length(test_aux_dwns))'/fsl);
    mesh(X, Y, crit_hist_viz, 'facealpha',0.5); hold on
    
    %Plot criterion maximal points
    [crit_maxs, max_inds] = max(crit_hist_viz,[],2);
    plot3(crit_nodes(max_inds), (1:length(test_aux_dwns))'/fsl, crit_maxs, 'r*');
    ylabel('Time, s'); xlabel('Effect $e$, normalized'); zlabel('$\log C_I$')
end

clear plot_criterion X Y

%%
%clear signals_file annotation_file fe_* clear EMGSIGNAL AUXSIGNAL current_path
%clear test_timeline test_aux* test_predicted* test_activity* i
%clear train_* test_*