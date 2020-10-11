%% mikexcohen.com
% phase-based connectivity

load sampleEEGdata

%% setup some stuff and junk

% names of the channels you want to compute connectivity between
channel1 = 'p1';
channel2 = 'pz';


% create complex Morlet wavelet
center_freq = 5;
time      = -1:1/EEG.srate:1;
wavelet   = exp(2*1i*pi*center_freq.*time) .* exp(-time.^2./(2*(4/(2*pi*center_freq))^2));
half_wavN = (length(time)-1)/2;

% FFT parameters
n_wavelet = length(time);
n_data    = EEG.pnts;
n_conv    = n_wavelet+n_data-1;

% FFT of wavelet
waveletX = fft(wavelet,n_conv);
waveletX = waveletX ./ max(waveletX);

% initialize output time-frequency data
phase_data = zeros(2,EEG.pnts);
real_data  = zeros(2,EEG.pnts);

% find channel indices
chan1idx = find(strcmpi(channel1,{EEG.chanlocs.labels}));
chan2idx = find(strcmpi(channel2,{EEG.chanlocs.labels}));


% analytic signal of channel 1
fft_data = fft(squeeze(EEG.data(chan1idx,:,1)),n_conv);
as = ifft(waveletX.*fft_data,n_conv);
as = as(half_wavN+1:end-half_wavN);

% collect real and phase data
phase_data(1,:) = angle(as);
real_data(1,:)  = real(as);

% analytic signal of channel 1
fft_data = fft(squeeze(EEG.data(chan2idx,:,1)),n_conv);
as = ifft(waveletX.*fft_data,n_conv);
as = as(half_wavN+1:end-half_wavN);

% collect real and phase data
phase_data(2,:) = angle(as);
real_data(2,:)  = real(as);

%% setup figure and define plot handles

% open and name figure
figure, set(gcf,'Name','Movie magic minimizes the magic.');

% draw the filtered signals
subplot(221)
filterplotH1 = plot(EEG.times(1),real_data(1,1),'b');
hold on
filterplotH2 = plot(EEG.times(1),real_data(2,1),'m');
set(gca,'xlim',[EEG.times(1) EEG.times(end)],'ylim',[min(real_data(:)) max(real_data(:))])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title([ 'Filtered signal at ' num2str(center_freq) ' Hz' ])

% draw the phase angle time series
subplot(222)
phaseanglesH1 = plot(EEG.times(1),phase_data(1,1),'b');
hold on
phaseanglesH2 = plot(EEG.times(1),phase_data(2,1),'m');
set(gca,'xlim',[EEG.times(1) EEG.times(end)],'ylim',[-pi pi]*1.1)
xlabel('Time (ms)')
ylabel('Phase angle (radian)')
title('Phase angle time series')

% draw phase angles in polar space
subplot(223)
polar2chanH1 = polar([zeros(1,1) phase_data(1,1)]',repmat([0 1],1,1)','b');
hold on
polar2chanH2 = polar([zeros(1,1) phase_data(2,1)]',repmat([0 1],1,1)','m');
title('Phase angles from two channels')

% draw phase angle differences in polar space
subplot(224)
polarAngleDiffH = polar([zeros(1,1) phase_data(2,1)-phase_data(1,1)]',repmat([0 1],1,1)','k');
title('Phase angle differences from two channels')

%% now update plots at each timestep

for ti=1:5:EEG.pnts
    
    % update filtered signals
    set(filterplotH1,'XData',EEG.times(1:ti),'YData',real_data(1,1:ti))
    set(filterplotH2,'XData',EEG.times(1:ti),'YData',real_data(2,1:ti))
    
    % update cartesian plot of phase angles
    set(phaseanglesH1,'XData',EEG.times(1:ti),'YData',phase_data(1,1:ti))
    set(phaseanglesH2,'XData',EEG.times(1:ti),'YData',phase_data(2,1:ti))
    
    subplot(223)
    cla
    polar([zeros(1,ti) phase_data(1,1:ti)]',repmat([0 1],1,ti)','b');
    hold on
    polar([zeros(1,ti) phase_data(2,1:ti)]',repmat([0 1],1,ti)','m');
    
    subplot(224)
    cla
    polar([zeros(1,ti) phase_data(2,1:ti)-phase_data(1,1:ti)]',repmat([0 1],1,ti)','k');
    
    drawnow
end

%% now quantify phase synchronization between the two channels

% phase angle differences
phase_angle_differences = phase_data(2,:)-phase_data(1,:);

% euler representation of angles
euler_phase_differences = exp(1i*phase_angle_differences);

% mean vector (in complex space)
mean_complex_vector = mean(euler_phase_differences);

% length of mean vector (this is the "M" from Me^ik, and is the measure of phase synchronization)
phase_synchronization = abs(mean_complex_vector);

disp([ 'Synchronization between ' channel1 ' and ' channel2 ' is ' num2str(phase_synchronization) '!' ])

% of course, this could all be done on one line:
phase_synchronization = abs(mean(exp(1i*(phase_data(2,:)-phase_data(1,:)))));

% notice that the order of subtraction is meaningless (see below), which means that this measure of synchronization is non-directional!
phase_synchronization_backwards = abs(mean(exp(1i*(phase_data(1,:)-phase_data(2,:)))));


% now plot mean vector
subplot(224)
hold on
h=polar([0 angle(mean_complex_vector)],[0 phase_synchronization]);
set(h,'linewidth',6,'color','g')

%% phase clustering is phase-invariant

figure(2), clf
subplot(221)
polar(repmat(phase_data(2,:)-phase_data(1,:),1,2)',repmat([0 1],1,EEG.pnts)','k');
title([ 'Phase synchronization: ' num2str(abs(mean(exp(1i*(diff(phase_data,1)))))) ])

new_phase_data = phase_data;
for i=2:4
    subplot(2,2,i)
    
    % add random phase offset
    new_phase_data(1,:) = new_phase_data(1,:)+rand*pi;
    
    % plot again
    polar(repmat(new_phase_data(2,:)-new_phase_data(1,:)+pi/2,1,2)',repmat([0 1],1,EEG.pnts)','k');
    title([ 'Phase synchronization: ' num2str(abs(mean(exp(1i*(diff(new_phase_data,1)))))) ])
end

%% 

%% ISPC over time vs. over trials

% FFT parameters
n_wavelet = length(time);
n_data    = EEG.pnts*EEG.trials;
n_conv    = n_wavelet+n_data-1;


% initialize output time-frequency data
phase_data = zeros(2,EEG.pnts,EEG.trials);


% FFT of wavelet (need to redo FFT because different n_conv)
waveletX = fft(wavelet,n_conv);
waveletX = waveletX ./ max(waveletX);


% analytic signal of channel 1
fft_data = fft(reshape(EEG.data(chan1idx,:,:),1,[]),n_conv);
as = ifft(waveletX.*fft_data,n_conv);
as = as(half_wavN+1:end-half_wavN);
as = reshape(as,EEG.pnts,EEG.trials);

% collect real and phase data
phase_data(1,:,:) = angle(as);

% analytic signal of channel 1
fft_data = fft(reshape(EEG.data(chan2idx,:,:),1,[]),n_conv);
as = ifft(waveletX.*fft_data,n_conv);
as = as(half_wavN+1:end-half_wavN);
as = reshape(as,EEG.pnts,EEG.trials);

% collect real and phase data
phase_data(2,:,:) = angle(as);


figure(3), clf

% ISPC over trials
subplot(211)
ispc_trials = abs(mean(exp(1i*diff(phase_data,[],1)
plot(EEG.times,ispc_trials)
set(gca,'xlim',[-200 1200])
xlabel('Time (ms)'), ylabel('ISPC')

% ISPC over time
subplot(212)
ispc_time = squeeze( abs(mean(exp(1i*diff(phase_data,[],1)
plot(1:EEG.trials,ispc_time)
xlabel('Trials'), ylabel('ISPC')

%% end.
