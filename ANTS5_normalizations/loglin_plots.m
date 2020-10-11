% mikexcohen@gmail.com

%% initial parameters

% vector of time points to save in post-analysis downsampling
times2save = -300:20:1200; % in ms

basetime = [-500 -200];

% frequency parameters
min_freq =  2;
max_freq = 50;
num_frex = 40;
% other wavelet parameters
range_cycles = [ 4 10 ];

% which channel to plot
channel2use = 'pz';

%% parameter conversion and other initializations

load sampleEEGdata.mat

% time vector converted to indices
times2saveidx = dsearchn(EEG.times',times2save');
basetimeidx   = dsearchn(EEG.times',basetime');

% frequencies vector
frex = logspace(log10(min_freq),log10(max_freq),num_frex);

% wavelet parameters
s = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex) ./ (2*pi*frex);
wavtime = -2:1/EEG.srate:2;
half_wave = (length(wavtime)-1)/2;


% FFT parameters
nWave = length(wavtime);
nData = EEG.pnts * EEG.trials; % This line is different from above!!
nConv = nWave + nData - 1;

% initialize output time-frequency data
tf = zeros(length(frex),length(times2save));

%% run convolution

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
    temppow  = mean( abs(as).^2 ,2);
    tf(fi,:) = 10*log10( temppow(times2saveidx) / mean(temppow(basetimeidx(1):basetimeidx(2))) );
end

%% plotting

clim = [-3 3];

figure(2)
subplot(121)
contourf(times2save,frex,tf,40,'linecolor','none')
set(gca,'clim',clim,'yscale','log')

% two different options for plotting
set(gca,'ytick',1:3:num_frex)
set(gca,'ytick',ceil(logspace(log10(1),log10(num_frex),8)))

title('Logarithmic frequency scaling')
xlabel('Time (ms)'), ylabel('Frequencies (Hz)')

subplot(122)
contourf(times2save,frex,tf,40,'linecolor','none')
set(gca,'clim',clim)
title('Linear frequency scaling')

%% imagesc vs. contourf

figure(2), clf
subplot(221)
contourf(times2save,frex,tf,40,'linecolor','none')
set(gca,'clim',clim,'xlim',[-200 1000],'yscale','log','ytick',ceil(logspace(log10(1),log10(num_frex),8)))
title('Logarithmic frequency scaling')
ylabel('Frequency (Hz)')

subplot(222)
contourf(times2save,frex,tf,40,'linecolor','none')
set(gca,'clim',clim,'xlim',[-200 1000])
title('Linear frequency scaling')

subplot(223)
imagesc(times2save,frex,tf)
set(gca,'clim',clim,'xlim',[-200 1000],'ydir','norm')
title('WRONG Y-AXIS LABELS!!!!')
ylabel('Frequency (Hz)'), xlabel('Time (ms)')

subplot(224)
imagesc(times2save,[],tf)
set(gca,'clim',clim,'xlim',[-200 1000],'ydir','norm')
set(gca,'ytick',round(linspace(1,num_frex,6)),'yticklabel',round(frex(round(linspace(1,num_frex,6)))*10)/10)
title('CORRECT Y-AXIS LABELS!!!!')
xlabel('Time (ms)')

%% end
