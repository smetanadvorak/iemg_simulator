function EMG = EMG_filtering_fir(EMG, low_freq, high_freq, notch_freq, fs)
%Filters one-channel EMG:
%EMG_filtering(EMG, lowest_freq, highest_freq, notch_freq (zero if no notch needed), sampling rate)

if low_freq

    hp_filt = designfilt('highpassfir',...
                     'StopbandFrequency', max(1,low_freq-fs/2/25),...
                     'PassbandFrequency', low_freq, ...
                     'StopbandAttenuation', 20, 'PassbandRipple', 0.5, ...
                     'SampleRate', fs, 'DesignMethod', 'equiripple');
    EMG = filtfilt(hp_filt, EMG);
    
end

if high_freq ~= inf

    hp_filt = designfilt('lowpassfir',...
                     'StopbandFrequency', max(1,high_freq+fs/2/25),...
                     'PassbandFrequency', high_freq, ...
                     'StopbandAttenuation', 20, 'PassbandRipple', 0.5, ...
                     'SampleRate', fs, 'DesignMethod', 'equiripple');
                 
    EMG = filtfilt(hp_filt, EMG);
end


for i = 1:numel(notch_freq)
        bs_filt = designfilt('bandstopfir','FilterOrder',fs/4, ...
            'CutoffFrequency1', notch_freq(i)-fs/500 ,'CutoffFrequency2', notch_freq(i)+fs/500, ...
            'SampleRate',fs);
        EMG = filtfilt(bs_filt, EMG);
end






