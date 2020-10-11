%% mikexcohen.com

%% the FFT-based hilbert transform

% generate random numbers
n = 21;
randomnumbers = randn(n,1);

% take FFT
f = fft(randomnumbers);
% create a copy that is multiplied by the complex operator
complexf = 1i*f;

% find indices of positive and negative frequencies
posF = 2:floor(n/2)+mod(n,2);
negF = ceil(n/2)+1+~mod(n,2):n;

% rotate Fourier coefficients
% (note 1: this works by computing the iAsin(2pft) component, i.e., the phase quadrature)
% (note 2: positive frequencies are rotated counter-clockwise; negative frequencies are rotated clockwise)
f(posF) = f(posF) + -1i*complexf(posF);
f(negF) = f(negF) +  1i*complexf(negF);
% The next two lines are an alternative and slightly faster method. 
% The book explains why this is equivalent to the previous two lines.
% f(posF) = f(posF)*2;
% f(negF) = f(negF)*0;

% take inverse FFT
hilbertx = ifft(f);

% compare with Matlab function hilbert
hilbertm = hilbert(randomnumbers);

% plot results
figure(1), clf
subplot(211)
plot(abs(hilbertm))
hold on
plot(abs(hilbertx),'ro')
set(gca,'xlim',[.5 n+.5])
legend({'Matlab Hilbert function';'"manual" Hilbert'})
title('magnitude of Hilbert transform')

subplot(212)
plot(angle(hilbertm))
hold on
plot(angle(hilbertx),'ro')
set(gca,'xlim',[.5 n+.5])
legend({'Matlab Hilbert function';'"manual" Hilbert'})
title('phase of Hilbert transform')

%% try with real data

load sampleEEGdata

% ERP, and its hilbert transform
erp  = squeeze(mean(EEG.data(48,:,:),3));
erpH = hilbert(erp);

figure(2), clf

% plot ERP and real part of Hilbert transformed ERP
subplot(311)
plot(EEG.times,erp), hold on
plot(EEG.times,real(erpH),'r')
legend({'ERP';'real(hilbert(erp))'})

% plot ERP and magnitude
subplot(312)
plot(EEG.times,erp), hold on
plot(EEG.times,abs(erpH),'r')
legend({'ERP';'abs(hilbert(erp))'})

% plot ERP and phase angle time series
subplot(313)
plot(EEG.times,erp), hold on
plot(EEG.times,angle(erpH),'r')
legend({'ERP';'angle(hilbert(erp))'})
xlabel('Time (ms)'), ylabel('Voltage or radians')



% plot as 3d line
figure(3), clf
plot3(EEG.times,real(erpH),imag(erpH))
xlabel('Time (ms)'), ylabel('Real part'), zlabel('Imaginary part')
axis tight
rotate3d

%% done.

