% mikexcohen@gmail.com

%% specify parameters and settings

load sampleEEGdata.mat

% select channels, time windows, and frequency ranges
channel1 = 'o1';
timewin1 = [ -250 0 ]; % in ms
freqwin1 = [ 8 12 ]; % in Hz

channel2 = 'fz';
timewin2 = [ 250 500 ]; % in ms
freqwin2 = [ 4 8 ]; % in Hz


% frequency parameters for time-frequency decomposition
min_freq =  2;
max_freq = 50;
num_frex = 40;
frex = linspace(min_freq,max_freq,num_frex);

% other wavelet parameters
range_cycles = [ 4 10 ];

%% convert connectivity parameters to indices

chan1idx = strcmpi({EEG.chanlocs.labels},channel1);
time1idx = dsearchn(EEG.times',timewin1');
freq1idx = dsearchn(frex',freqwin1');

chan2idx = strcmpi({EEG.chanlocs.labels},channel2);
time2idx = dsearchn(EEG.times',timewin2');
freq2idx = dsearchn(frex',freqwin2');

%% setup wavelet and convolution parameters

s = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex) ./ (2*pi*frex);
wavtime = -2:1/EEG.srate:2;
half_wave = (length(wavtime)-1)/2;


% convolution parameters
nWave = length(wavtime);
nData = EEG.pnts * EEG.trials;
nConv = nWave + nData - 1;

% initialize output time-frequency data
tfall = zeros(2,length(frex),EEG.pnts,EEG.trials);
data4corr = zeros(2,EEG.trials);

%% run convolution

for chani=1:2
    
    % now compute the FFT of all trials concatenated
    if chani==1
        alldata = reshape( EEG.data(chan1idx,:,:) ,1,[]);
    else 
        alldata = reshape( EEG.data(chan2idx,:,:) ,1,[]);
    end
    
    % this line does the same as the previous if-else statement
    eval([ 'alldata = reshape( EEG.data(chan' num2str(chani) 'idx,:,:) ,1,[]);' ]);
    
    % Fourier spectrum of data
    dataX = fft( alldata ,nConv );
    
    % loop over frequencies
    for fi=1:length(frex)
        
        % create wavelet and get its FFT
        % (does this need to be computed here inside this double-loop?)
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv);
        waveletX = waveletX ./ max(waveletX);
        
        % now run convolution in one step
        as = ifft(waveletX .* dataX);
        as = as(half_wave+1:end-half_wave);
        
        % and reshape back to time X trials
        as = reshape( as, EEG.pnts, EEG.trials );
        
        % compute power and save for all trials
        tfall(chani,fi,:,:) = abs(as).^2;
    end
end

%% check that the data look OK

figure(1), clf

for chani=1:2
    subplot(2,1,chani)
    contourf(EEG.times,frex,squeeze(mean( tfall(chani,:,:,:) ,4)),40,'linecolor','none')
    set(gca,'clim',[0 5],'xlim',[-200 1300])
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    eval([ 'title([ ''Channel '' channel' num2str(chani) ' ])' ])
end

%% and now extract data for connectivity

data4corr(1,:) = squeeze(mean(mean( tfall(1,freq1idx(1):freq1idx(2),time1idx(1):time1idx(2),:) ,2),3));
data4corr(2,:) = squeeze(mean(mean( tfall(2,freq2idx(1):freq2idx(2),time2idx(1):time2idx(2),:) ,2),3));

figure(2), clf
plot(data4corr(1,:),data4corr(2,:),'o')
xlabel([ 'Power: ' channel1 ', ' num2str(freqwin1(1)) '-' num2str(freqwin1(2)) 'Hz, ' num2str(timewin1(1)) '-' num2str(timewin1(2)) 'ms' ])
ylabel([ 'Power: ' channel2 ', ' num2str(freqwin2(1)) '-' num2str(freqwin2(2)) 'Hz, ' num2str(timewin2(1)) '-' num2str(timewin2(2)) 'ms' ])

[r,p] = corr(data4corr','type','spearman');

title([ 'Correlation R=' num2str(r(1,2)) ', p=' num2str(p(1,2)) ])

%% end
