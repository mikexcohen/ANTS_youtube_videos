% mikexcohen@gmail.com

%% show power spectrum of EEG data

load sampleEEGdata.mat

% pick a channel
chan2plot = 'tp8';


% compute Fourier spectrum
eegX = fft( squeeze(EEG.data(strcmpi(chan2plot,{EEG.chanlocs.labels}),:,:)) ,[],1)/EEG.pnts;
ampl = 2*abs(eegX);


% Compute frequencies in Hz. Note that this is a different way of compute
% frequencies, because frequencies are specified up to the sampling rate,
% not the Nyquist frequency! It's worth knowing this trick, but best to
% avoid it most of the time.
hz = linspace(0,EEG.srate,EEG.pnts);

% and plot
figure(1), clf
plot(hz,ampl)
hold on
plot(hz,mean(ampl,2),'k','linew',4)

set(gca,'xlim',[0 55])
xlabel('Frequency (Hz)'), ylabel('Amplitude (\muV)')
title([ 'Single-trial and trial-average power spectrum from channel ' chan2plot ])

% 
% set(gca,'yscale','log')
% set(gca,'xscale','log')

%% end
