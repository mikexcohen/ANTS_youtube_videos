% mikexcohen@gmail.com

%% TF power decomposition

load sampleEEGdata.mat

% vector of time points to save in post-analysis downsampling
times2save = -300:20:1200; % in ms


% time vector converted to indices
times2saveidx = dsearchn(EEG.times',times2save');


% frequency parameters
min_freq =  2;
max_freq = 50;
num_frex = 40;
frex = linspace(min_freq,max_freq,num_frex);

% which channel to plot
channel2use = 'o1';

% other wavelet parameters
range_cycles = [ 4 10 ];

s = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex) ./ (2*pi*frex);
wavtime = -2:1/EEG.srate:2;
half_wave = (length(wavtime)-1)/2;


% FFT parameters
nWave = length(wavtime);
nData = EEG.pnts * EEG.trials; % This line is different from above!!
nConv = nWave + nData - 1;

% initialize output time-frequency data
tf = zeros(length(frex),EEG.pnts);

% now compute the FFT of all trials concatenated
alldata = reshape( EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:) ,1,[]);
dataX   = fft( alldata ,nConv );


% loop over frequencies
for fi=1:length(frex)
    
    % create wavelet and get its FFT
    % the wavelet doesn't change on each trial...
    wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
    waveletX = fft(wavelet,nConv);
    waveletX = waveletX ./ max(waveletX);
    
    % now run convolution in one step
    as = ifft(waveletX .* dataX);
    as = as(half_wave+1:end-half_wave);
    
    % and reshape back to time X trials
    as = reshape( as, EEG.pnts, EEG.trials );
    
    % compute power and average over trials
    tf(fi,:) = mean( abs(as).^2 ,2);
end

% plot results
figure(1), clf
subplot(211)
contourf(EEG.times,frex,tf,40,'linecolor','none')
set(gca,'clim',[0 5],'ydir','normal','xlim',[-300 1000])
title('Full temporal resolution')


subplot(212)
contourf(times2save,frex,tf(:,times2saveidx),40,'linecolor','none')
set(gca,'clim',[0 5],'ydir','normal','xlim',[-300 1000])
title('Temporally down-sampled after convolution')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')


% now show lines
figure(2), clf
% pick frequencies to isolate
freqs2plot = dsearchn(frex',[5 30]');
plot(EEG.times,tf(freqs2plot,:),'.-')
hold on
plot(times2save,tf(freqs2plot,times2saveidx),'o-')
xlabel('Time (ms)'), ylabel('Power (\muv^2)')

legend({[ num2str(round(frex(freqs2plot(1)))) ' Hz full' ];[ num2str(round(frex(freqs2plot(2)))) ' Hz full' ]; ...
        [ num2str(round(frex(freqs2plot(1)))) ' Hz downsampled' ];[ num2str(round(frex(freqs2plot(2)))) ' Hz downsampled' ]})

%% how to downsample during real analyses

% initialize output time-frequency data
tfds = zeros(length(frex),length(times2save));

% loop over frequencies
for fi=1:length(frex)
    
    % create wavelet and get its FFT
    % the wavelet doesn't change on each trial...
    wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
    waveletX = fft(wavelet,nConv);
    waveletX = waveletX ./ max(waveletX);
    
    % now run convolution in one step
    as = ifft(waveletX .* dataX);
    as = as(half_wave+1:end-half_wave);
    
    % and reshape back to time X trials
    as = reshape( as, EEG.pnts, EEG.trials );
    
    % compute power and average over trials
    % here is where we take only the downsampled times
    tfds(fi,:) = mean( abs(as(times2saveidx,:)).^2 ,2);
end

whos tf*

%% end
