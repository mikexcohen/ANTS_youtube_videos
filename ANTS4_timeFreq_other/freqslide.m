%% mikexcohen.com

%% simulate FM signal

% signal properties
srate = 1000;
time  = 0:1/srate:5;

% time series of instantaneous frequencies, which is created by
% interpolating a few random datapoints
freqTS = abs( interp1(linspace(time(1),time(end),10),10*rand(1,10),time,'spline') );

% this is the construction of a frequency-varying sine wave
centfreq = mean(freqTS);
k        = (centfreq/srate)*2*pi/centfreq;
signal   = sin(2*pi.*centfreq.*time + k*cumsum(freqTS-centfreq));

% hilbert transform of signal
hilsig = hilbert(signal);


% and plot!
figure(1), clf
subplot(311)
plot(time,signal)
ylabel('Amplitude')
set(gca,'ylim',[-1.1 1.1])

subplot(312)
plot(time,angle(hilsig))
ylabel('Phase (rad.)')

subplot(313)
freqslide = srate*diff(unwrap(angle(hilsig)))/(2*pi);
plot(time,freqslide)
hold on
plot(time,freqTS,'r')
set(gca,'ylim',[0 15])
legend({'Estimated';'Simulated'})

%% And now with real data and a median filter

load sampleEEGdata

data2use = squeeze(EEG.data(28,:,10));
freq2use = 10; % hz


% define convolution parameters
wavt  = -2:1/EEG.srate:2; % time vector for wavelet
nData = EEG.pnts;
nKern = length(wavt);
nConv = nData+nKern-1;
hwave = floor((length(wavt)-1)/2);

s = 8 / (2*pi*freq2use); % for gaussian
cmwX = fft( exp(1i*2*pi*freq2use*wavt) .* exp( -(wavt.^2)/(2*s^2) ) ,nConv);
cmwX = cmwX ./ max(cmwX);

as = ifft( fft(data2use,nConv) .* cmwX );
as = as(hwave+1:end-hwave);


figure(2), clf
subplot(311)
plot(EEG.times,data2use);
ylabel('Amplitude (\muV)')

subplot(312)
plot(EEG.times,angle(as))
ylabel('Phase angles (rad.)')

subplot(313)
freqslide = EEG.srate*diff(unwrap(angle(as)))/(2*pi);
plot(EEG.times(1:end-1),freqslide)
ylabel('Frequency (Hz)')
set(gca,'ylim',[freq2use*.75 freq2use*1.5])


%% now apply median filter

n_order = 10;
orders  = linspace(10,400,n_order)/2; 
orders  = round( orders/(1000/EEG.srate) );


phasedmed = zeros(length(orders),EEG.pnts);

for oi=1:n_order
    for ti=1:EEG.pnts
        
        %% use compiled fast_median if available
%         phasedmed(oi,ti) = fast_median(freqslide(max(ti-orders(oi),1):min(ti+orders(oi),EEG.pnts-1))');
        
        %% use 'manual median' otherwise
        temp = sort(freqslide(max(ti-orders(oi),1):min(ti+orders(oi),EEG.pnts-1)));
        phasedmed(oi,ti) = temp(floor(numel(temp)/2)+1);
        
    end
end

% the final step is to take the mean of medians
freqslideFilt = mean(phasedmed,1);

subplot(313)
hold on
plot(EEG.times,freqslideFilt,'r')

%%
