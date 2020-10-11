%% mikexcohen.com

%% specify baseline periods for dB-conversion

load sampleEEGdata

baseline_windows = [ -500 -200;
                     -100    0;
                        0  300;
                     -800    0;
                   ];

               
% convert baseline time into indices
baseidx = reshape( dsearchn(EEG.times',baseline_windows(:)), [],2);

%% setup wavelet parameters

% frequency parameters
min_freq =  2;
max_freq = 30;
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
nData = EEG.pnts * EEG.trials;
nConv = nWave + nData - 1;


% now compute the FFT of all trials concatenated
alldata = reshape( EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:) ,1,[]);
dataX   = fft( alldata ,nConv );


% initialize output time-frequency data
tf = zeros(size(baseidx,1),length(frex),EEG.pnts);

%% now perform convolution


% loop over frequencies
for fi=1:length(frex)
    
    % create wavelet and get its FFT
    wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
    waveletX = fft(wavelet,nConv);
    waveletX = waveletX ./ max(waveletX);
    
    % now run convolution in one step
    as = ifft(waveletX .* dataX);
    as = as(half_wave+1:end-half_wave);
    
    % and reshape back to time X trials
    as = reshape( as, EEG.pnts, EEG.trials );
    
    % compute power and average over trials
    tf(4,fi,:) = mean( abs(as).^2 ,2);
end

%% db conversion and plot results

% define color limits
climdb  = [-3 3];
climpct = [-90 90];

% create new matrix for percent change
tfpct = zeros(size(tf));

for basei=1:size(tf,1)
    
    activity = tf(4,:,:);
    baseline = mean( tf(4,:,baseidx(basei,1):baseidx(basei,2)) ,3);
    
    % decibel
    tf(basei,:,:) = 10*log10( bsxfun(@rdivide, activity, baseline) );
    
    % percent change
    tfpct(basei,:,:) = 100 * bsxfun(@rdivide, bsxfun(@minus,activity,baseline), baseline);
end


% plot results
for basei=1:size(baseline_windows,1)
    
    % first plot dB
    figure(1), subplot(2,2,basei)
    
    contourf(EEG.times,frex,squeeze(tf(basei,:,:)),40,'linecolor','none')
    set(gca,'clim',climdb,'ydir','normal','xlim',[-300 1000])
    title([ 'DB baseline of ' num2str(baseline_windows(basei,1)) ' to ' num2str(baseline_windows(basei,2)) ' ms' ])
    
    % now plot percent change
    figure(2), subplot(2,2,basei)
    
    contourf(EEG.times,frex,squeeze(tfpct(basei,:,:)),40,'linecolor','none')
    set(gca,'clim',climpct,'ydir','normal','xlim',[-300 1000])
    title([ 'PCT baseline of ' num2str(baseline_windows(basei,1)) ' to ' num2str(baseline_windows(basei,2)) ' ms' ])
end

xlabel('Time (ms)'), ylabel('Frequency (Hz)')

%%
