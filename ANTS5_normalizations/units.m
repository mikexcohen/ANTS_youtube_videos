% mikexcohen.com

%% units in FFT amplitude

% define parameters
srate = 1000;
time  = 0:1/srate:10;
nPnts = length(time);
freq  = 14; % Hz
hz    = linspace(0,srate/2,floor(nPnts/2)+1);

% create signal
sinewave = 2*sin(2*pi*freq*time);

% step 1: FFT and divide by N
sinewaveX = fft(sinewave,nPnts)/nPnts; % divide by number of points

% step 2: multiply positive-frequency coefficients by 2
posfrex  = 2:floor(nPnts/2);
sineampl = abs(sinewaveX);
sineampl(posfrex) = sineampl(posfrex)*2;

figure(1), clf
subplot(211)
plot(time,sinewave)
set(gca,'xlim',[0 1])
xlabel('Time (s)'), ylabel('Amplitude (a.u.)')

subplot(212)
bar(hz,sineampl(1:length(hz)))
set(gca,'xlim',[freq-3 freq+3])
xlabel('Frequency (Hz)'), ylabel('Amplitude (a.u.)')

%% units in FFT power

% define parameters
srate = 1000;
time  = 0:1/srate:10;
nPnts = length(time);
freq  = 14; % Hz
hz    = linspace(0,srate/2,floor(nPnts/2)+1);

% create signal
sinewave = 2*sin(2*pi*freq*time);

% step 1: FFT and divide by N
sinewaveX = fft(sinewave)/nPnts; % divide by 

% step 2: multiply positive-frequency coefficients by 2
posfrex  = 2:floor(nPnts/2);
sineampl = abs(sinewaveX);
sineampl(posfrex) = sineampl(posfrex)*2;

% step 3: square amplitudes
sineampl = sineampl.^2;


figure(2), clf
subplot(211)
plot(time,sinewave)
set(gca,'xlim',[0 1])
xlabel('Time (s)'), ylabel('Amplitude (a.u.)')

subplot(212)
bar(hz,sineampl(1:length(hz)))
set(gca,'xlim',[freq-3 freq+3])
xlabel('Power (Hz)'), ylabel('Amplitude (a.u.)')

%% setting up wavelet convolution in preparation for extracting correct units

% create signal
srate  = 1000;
time   = 0:1/srate:10;
nPnts  = length(time);
frex   = [25 60];
signal = 2*sin(2*pi.*linspace(frex(1),frex(2)*mean(frex)/frex(2),nPnts).*time);


% specify wavelet parameters
wfreq  = 40;               % in Hz
nCyc   = 20/(2*pi*wfreq);  % numerator is number of cycles; denominator is a scaling factor

% create wavelet
wavet   = -1:1/srate:1;
srate   = 1000;
halfwav = (length(wavet)-1)/2;
cmw     = exp(1i*2*pi*wfreq*wavet) .* exp( -.5*(wavet/nCyc).^2 );

% convolution parameters
nData = nPnts;
nKern = length(wavet);
nConv = nData+nKern-1;

% FFT of wavelet and signal
cmwX  = fft(cmw,nConv);
sigX  = fft(signal,nConv);

% step 1: max-value normalize the wavelet in the frequency domain
cmwX  = cmwX./max(cmwX);

% convolution
asX   = ifft( cmwX.*sigX );

% frequencies in Hz (valid only up to Nyquist)
hz    = -(srate/2-1/srate):srate/nConv:srate/2;

% plotting
figure(3), clf
subplot(311), plot(time,signal)
set(gca,'xlim',[2 5])
title('Signal in the time domain')
xlabel('Time (sec)'), ylabel('Amplitude (a.u.)')

subplot(312)
plot(hz,fftshift(2*abs(sigX/nPnts)))
title('Signal in the frequency domain')
xlabel('Frequency (Hz)'), ylabel('Amplitude (normalized)')

subplot(313)
plot(hz,fftshift(2*abs(cmwX)))
title('Wavelet in the frequency domain')
xlabel('Frequency (Hz)'), ylabel('Amplitude (normalized)')

%% now get units back

% trim convolution wings
asX = asX(halfwav+1:end-halfwav);

tf_amp = 2*abs(asX);
tf_pow = tf_amp.^2;


figure(4), clf
subplot(211)
plot(time,signal), hold on
plot(time,tf_amp,'linew',2)
xlabel('Time (s)'), ylabel('Amplitude (a.u.)')

subplot(212)
plot(time,signal), hold on
plot(time,tf_pow,'linew',2)
xlabel('Time (s)'), ylabel('Power (a.u.)')

%% end.
