%% mikexcohen.com

%% short-window FFT

load sampleEEGdata.mat

chan2plot  = 'o1';

frex       = logspace(log10(10),log10(EEG.srate/5),20);
times2save = -300:25:800;
basetime   = [-300 -100];
timewin    = 300; % in ms


% convert time points to indices
times2saveidx = dsearchn(EEG.times',times2save'); 
% convert time window to points
timewinpnts   = round(timewin/(1000/EEG.srate));

% find baselinetimepoints
baseidx = dsearchn(times2save',basetime');

% define frequencies for FFT
hz = linspace(0,EEG.srate/2,timewinpnts/2+1);

% hanning window for tapering
hannwin = .5 - .5*cos(2*pi.*linspace(0,1,timewinpnts))';

% find logical channel index
chanidx = strcmpi(chan2plot,{EEG.chanlocs.labels});

% initialize output matrix
shortFFT_tf = zeros(length(frex),length(times2save));


% loop through time bins
for ti=1:length(times2saveidx)
    
    % window and taper data, and get power spectrum
    data = bsxfun(@times, squeeze(EEG.data(chanidx,times2saveidx(ti)-floor(timewinpnts/2)+1:times2saveidx(ti)+ceil(timewinpnts/2),:)), hannwin);
    % uncomment the next line to use non-tapered data
    %data = squeeze(EEG.data(chanidx,times2saveidx(ti)-floor(timewinpnts/2)+1:times2saveidx(ti)+ceil(timewinpnts/2),:));
    
    y    = fft(data,timewinpnts)/timewinpnts;
    pow  = mean(abs(y).^2,2);
    
    % finally, get power from closest frequency
    closestfreq = dsearchn(hz',frex');
    shortFFT_tf(:,ti) = pow(closestfreq);    
end % end time loop

% db-correct
db_shortFFT_tf = 10*log10( bsxfun(@rdivide,shortFFT_tf,mean(shortFFT_tf(:,baseidx(1):baseidx(2)),2)) );

% plot!
figure(1), clf
contourf(times2save,frex,db_shortFFT_tf,40,'linecolor','none')
set(gca,'clim',[-2 2])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Power via short-window FFT (window=' num2str(timewin) ') from channel ' chan2plot ])

%%
