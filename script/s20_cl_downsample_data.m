
if fs == 10240
    fsl = 1024;
else
    fsl = 1000;
end


%% Define downsamplers
downsampler_aux = Downsampler(fs/fsl);
downsampler_spk = Downsampler(fs/fsl, size(train_spk,2));
downsampler_act = Downsampler(fs/fsl, size(train_spk,2));

actwin = round(250/1000 * fs);

%% Training data
train_aux_dwns = zeros(size(train_aux));
train_act_dwns = zeros(size(train_spk));
train_spk_dwns  = zeros(size(train_spk));
train_act = spikes2activity(train_spk, actwin, 1);

for t = 1:size(train_spk,1)
    train_aux_dwns(t) = downsampler_aux.update(train_aux(t));
    train_act_dwns(t,:) = downsampler_act.update(train_act(t,:));
    train_spk_dwns(t,:) = downsampler_spk.update(train_spk(t,:));
end

train_aux_dwns(isnan(train_aux_dwns(:,1)),:) = [];
train_act_dwns(isnan(train_act_dwns(:,1)),:) = [];
train_spk_dwns(isnan(train_spk_dwns(:,1)),:) = [];

train_act_dwns = double(train_act_dwns > 0);
train_spk_dwns = double(train_spk_dwns > 0);


%% Testing data
downsampler_aux.reset();
downsampler_spk.reset();
downsampler_act.reset();

test_aux_dwns = zeros(size(test_aux));
test_act_dwns = zeros(size(test_spk));
test_spk_dwns  = zeros(size(test_spk));
test_act = spikes2activity(test_spk, actwin, 1);

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