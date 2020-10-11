%% mikexcohen.com

%% band-pass filtering

load sampleEEGdata

% specify Nyquist freuqency
nyquist = EEG.srate/2;

% filter frequency band
filtbound = [4 10]; % Hz

% transition width
trans_width = 0.2; % fraction of 1, thus 20%

% filter order
filt_order = round(3*(EEG.srate/filtbound(1)));

% frequency vector (as fraction of Nyquist
ffrequencies  = [ 0 (1-trans_width)*filtbound(1) filtbound (1+trans_width)*filtbound(2) nyquist ]/nyquist;

% shape of filter (must be the same number of elements as frequency vector
idealresponse = [ 0 0 1 1 0 0 ];

% get filter weights
filterweights = firls(filt_order,ffrequencies,idealresponse);

% plot for visual inspection
figure(1), clf
subplot(211)
plot(ffrequencies*nyquist,idealresponse,'k--o','markerface','m')
set(gca,'ylim',[-.1 1.1],'xlim',[-2 nyquist+2])
xlabel('Frequencies (Hz)'), ylabel('Response amplitude')

subplot(212)
plot((0:filt_order)*(1000/EEG.srate),filterweights)
xlabel('Time (ms)'), ylabel('Amplitude')

% apply filter to data
filtered_data = zeros(EEG.nbchan,EEG.pnts);
for chani=1:EEG.nbchan
    filtered_data(chani,:) = filtfilt(filterweights,1,double(EEG.data(chani,:,1)));
end

figure(2), clf
plot(EEG.times,squeeze(EEG.data(47,:,1)))
hold on
plot(EEG.times,squeeze(filtered_data(chani,:)),'r','linew',2)
xlabel('Time (ms)'), ylabel('Voltage (\muV)')
legend({'raw data';'filtered'})

%% compare three transition widths

nyquist    = EEG.srate/2;
filtbond   = [ 7 12 ];
t_widths   = [.1 .15 .2];
filt_order = round(3*(EEG.srate/filtbond(1)));

idealresponse = [ 0 0 1 1 0 0 ];

ffrequencies  = zeros(3,6);
filterweights = zeros(3,filt_order+1);

% frequency vector (as fraction of Nyquist)
for i=1:3
    ffrequencies(i,:)  = [ 0 (1-t_widths(i))*filtbond(1) filtbond (1+t_widths(i))*filtbond(2) nyquist ]/nyquist;
    filterweights(i,:) = firls(filt_order,ffrequencies(i,:),idealresponse);
end

% plot
figure(4), clf
subplot(211)
plot((1:filt_order+1)*(1000/EEG.srate),filterweights')
xlabel('time (ms)')
title('Time-domain filter kernel')

filterFreqDomain = abs(fft(filterweights,[],2));
frequenciesHz    = linspace(0,nyquist,floor(filt_order/2)+1);
subplot(212)
plot(frequenciesHz,filterFreqDomain(:,1:length(frequenciesHz)))
set(gca,'xlim',[0 60],'ylim',[-.1 1.2])
xlabel('Frequencies (Hz)')
title('Frequency-domain filter kernel')
legend({'filter 10%','filter 15%','filter 20%'})

%% compute and plot power

chan4filt = strcmpi('o1',{EEG.chanlocs.labels});
baseidx   = dsearchn(EEG.times',[-400 -100]');

pow = zeros(3,EEG.pnts);

for i=1:3
    filtered_data = reshape(filtfilt(filterweights(i,:),1,double(reshape(EEG.data(chan4filt,:,:),1,[]))),EEG.pnts,EEG.trials);
    
    temppow  = mean(abs(hilbert(filtered_data)).^2,2);
    pow(i,:) = 10*log10( temppow./mean(temppow(baseidx(1):baseidx(2))) );
end

figure(5), clf
plot(EEG.times,pow)
xlabel('Time (ms)'), ylabel('power (dB)')
legend({'filter 10%','filter 15%','filter 20%'})
set(gca,'xlim',[-300 1200])

%%