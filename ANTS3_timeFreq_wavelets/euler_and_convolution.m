% mikexcohen.com


%% load in EEG data and pick and electrode, trial, etc.

load sampleEEGdata.mat

chan2use  = 'fz';
trial2use = 34;

% extract this bit of data for convenience in the rest of this script
data = squeeze(EEG.data( strcmpi(chan2use,{EEG.chanlocs.labels}), :,trial2use));

%% create a Morlet wavelet

srate = EEG.srate; % in hz
time  = -2:1/srate:2; % best practice is to have time=0 at the center of the wavelet
frex  = 6.5; % frequency of wavelet, in Hz

% create complex sine wave
sine_wave = exp( 1i*2*pi*frex.*time );

% create Gaussian window
s = 7 / (2*pi*frex); % this is the standard deviation of the gaussian
gaus_win  = exp( (-time.^2) ./ (2*s^2) );

% now create Morlet wavelet
cmw  = sine_wave .* gaus_win;

%% define convolution parameters

nData = length(data);
nKern = length(cmw);
nConv = nData + nKern - 1;
half_wav = floor( length(cmw)/2 )+1;

%% FFTs

% note that the "N" parameter is the length of convolution, NOT the length
% of the original signals! Super-important!


% FFT of wavelet, and amplitude-normalize in the frequency domain
cmwX = fft(cmw,nConv);
cmwX = cmwX ./ max(cmwX);

% FFT of data
dataX = fft(data,nConv);

% now for convolution...
as = ifft( dataX .* cmwX );

% cut 1/2 of the length of the wavelet from the beginning and from the end
as = as(half_wav-1:end-half_wav);

%% different ways of visualizing the outputs of convolution

figure(1), clf
plot3(EEG.times,real(as),imag(as),'k')
xlabel('Time (ms)'), ylabel('real part'), zlabel('imaginary part')
rotate3d
set(gca,'xlim',[-300 1200])

figure(2), clf
plot3(EEG.times,abs(as),angle(as),'k')
xlabel('Time (ms)'), ylabel('Amplitude'), zlabel('Phase')
rotate3d
set(gca,'xlim',[-300 1200])



figure(3), clf
% plot the filtered signal (projection onto real axis)
subplot(311)
plot(EEG.times,real(as))
xlabel('Time (ms)'), ylabel('Amplitude (\muV)')
set(gca,'xlim',[-200 1300])


% plot power (squared magnitude from origin to dot-product location in
% complex space)
subplot(312)
plot(EEG.times,abs(as).^2)
xlabel('Time (ms)'), ylabel('Power \muV^2')
set(gca,'xlim',[-200 1300])


% plot phase (angle of vector to dot-product, relative to positive real
% axis)
subplot(313)
plot(EEG.times,angle(as))
xlabel('Time (ms)'), ylabel('Phase (rad.)')
set(gca,'xlim',[-200 1300])

%% final view: in a polar plot

figure(4), clf

% setup the time course plot
subplot(212)
h = plot(EEG.times,abs(as),'k','linew',2);
set(gca,'xlim',EEG.times([1 end]),'ylim',[min(abs(as)) max(abs(as)) ])
xlabel('Time (ms)'), ylabel('Amplitude (\muV)')


for ti=1:2:EEG.pnts
    
    % draw complex values in polar space
    subplot(211)
    polar(0,max(abs(as))), hold on
    polar(angle(as(max(1,ti-100):ti)),abs(as(max(1,ti-100):ti)),'k')
    text(-.75,0,[ num2str(round(EEG.times(ti))) ' ms' ]), hold off
    
    % now show in 'linear' plot
    set(h,'XData',EEG.times(max(1,ti-100):ti),'YData',abs(as(max(1,ti-100):ti)))
    
    drawnow
    pause(.1)
end

%%
