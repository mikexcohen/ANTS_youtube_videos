%% mikexcohen.com

%% show tapers

load sampleEEGdata

timewin    = 400; % in ms
timewinidx = round(timewin/(1000/EEG.srate));
tapers     = dpss(timewinidx,5); % this line will crash without matlab signal processing toolbox

% plot tapers
figure(1), clf
for i=1:5
    subplot(5,1,i)
    plot( 1000*(0:timewinidx-1)/EEG.srate, tapers(:,i))
end
xlabel('Time (ms)')

%% Compute multitaper

timewin    = 300; % in ms
timewinidx = round(timewin/(1000/EEG.srate));
tapers     = dpss(timewinidx,3); % this line will crash without matlab signal processing toolbox

channel2plot = 'o1';
times2save   = -300:25:800;
basetime     = [-300 -100];

% convert time points to indices
times2saveidx = dsearchn(EEG.times',times2save'); 

% find baseline time point range
baseidx = dsearchn(times2save',basetime');

% define frequencies for FFT
hz = linspace(0,EEG.srate/2,timewinidx/2+1);

% find logical channel index
chanidx = strcmpi(channel2plot,{EEG.chanlocs.labels});

% initialize output matrix
tf = zeros(floor(timewinidx/2)+1,length(times2save));

% loop through time bins
for ti=1:length(times2saveidx)
    
    % initialize power vector (over tapers)
    taperpow = zeros(floor(timewinidx/2)+1,1);
    
    % loop through tapers
    for tapi = 1:size(tapers,2)-1
        
        % get data from this time window and taper
        tempEEG  = squeeze(EEG.data(chanidx,times2saveidx(ti)-floor(timewinidx/2)+1:times2saveidx(ti)+ceil(timewinidx/2),:));
        data     = bsxfun(@times,tempEEG,tapers(:,tapi));
        
        % compute FFT and extract power
        pow      = fft(data)/timewinidx;
        pow      = pow(1:floor(timewinidx/2)+1,:);
        taperpow = taperpow + mean(abs(pow).^2,2);
    end
    
    % divide by N tapers for average
    tf(:,ti) = taperpow/tapi;
end

% db-correct
tf = 10*log10( bsxfun(@rdivide,tf,mean(tf(:,baseidx(1):baseidx(2)),2)) );

% plot!
figure(2), clf
contourf(times2save,hz,tf,40,'linecolor','none')
set(gca,'clim',[-2 2],'ylim',[10 EEG.srate/5])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Power via multitaper from channel ' channel2plot ])

%% end.
