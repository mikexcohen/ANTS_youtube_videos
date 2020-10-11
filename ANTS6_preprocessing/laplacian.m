% mikexcohen.com

load sampleEEGdata

%% Simulated data

% which channels to use as hotspots
chan1 = 'pz';
chan2 = 'c4';
chan3 = 'c3';

% get XYZ coordinates in convenient variables
X = [EEG.chanlocs.X];
Y = [EEG.chanlocs.Y];
Z = [EEG.chanlocs.Z];

% initialize distance matrices
eucdist1 = zeros(1,64);
eucdist2 = zeros(1,64);
eucdist3 = zeros(1,64);

% convert electrode names to indices
chan1idx = strcmpi(chan1,{EEG.chanlocs.labels});
chan2idx = strcmpi(chan2,{EEG.chanlocs.labels});
chan3idx = strcmpi(chan3,{EEG.chanlocs.labels});

% compute inter-electrode distances, seeded from the three specified electrodes
for chani=1:EEG.nbchan
    eucdist1(chani) = sqrt( (EEG.chanlocs(chani).X-EEG.chanlocs(chan1idx).X)^2 + (EEG.chanlocs(chani).Y-EEG.chanlocs(chan1idx).Y)^2 + (EEG.chanlocs(chani).Z-EEG.chanlocs(chan1idx).Z)^2 );
    eucdist2(chani) = sqrt( (EEG.chanlocs(chani).X-EEG.chanlocs(chan2idx).X)^2 + (EEG.chanlocs(chani).Y-EEG.chanlocs(chan2idx).Y)^2 + (EEG.chanlocs(chani).Z-EEG.chanlocs(chan2idx).Z)^2 );
    eucdist3(chani) = sqrt( (EEG.chanlocs(chani).X-EEG.chanlocs(chan3idx).X)^2 + (EEG.chanlocs(chani).Y-EEG.chanlocs(chan3idx).Y)^2 + (EEG.chanlocs(chani).Z-EEG.chanlocs(chan3idx).Z)^2 );
end

% The code below defines the "activations," which are really 
% just a Gaussian function of distance away from each electrode.
% You can try changing widths to see how it affects the scalp maps.
lo_width = 95;
hi_width = 50;

lo_spatfreq  = 2*exp(- (eucdist1.^2)/(2*lo_width^2) ); 
hi_spatfreq  =   exp(- (eucdist2.^2)/(2*hi_width^2) )  +  exp(- (eucdist3.^2)/(2*hi_width^2) );


% and now compute the laplacian of the sum of the aforecreated data
surf_lap_all = laplacian_perrinX(hi_spatfreq+lo_spatfreq,X,Y,Z);



% show plots
figure(1), clf
subplot(221)
topoplot(lo_spatfreq,EEG.chanlocs,'plotrad',.53);
title('Low spatial frequency feature')

subplot(222)
topoplot(hi_spatfreq,EEG.chanlocs,'plotrad',.53);
title('High spatial frequency features')

subplot(223)
topoplot(lo_spatfreq+hi_spatfreq,EEG.chanlocs,'plotrad',.53);
title('Low+high features')

subplot(224)
topoplot(surf_lap_all,EEG.chanlocs,'plotrad',.53);
title('Laplacian of low+high features')

%% now for real data

times2plot = 0:200:800;

figure(2), clf
for i=1:length(times2plot)
    
    % find time index
    [~,timeidx] = min(abs(EEG.times-times2plot(i)));
    
    % average data from this time point
    tempdata = squeeze(mean(EEG.data(:,timeidx,:),3));
    
    
    % plot voltage map (spatially unfiltered)
    subplot(2,length(times2plot),i)
    topoplot(tempdata,EEG.chanlocs,'plotrad',.53,'maplimits',[-10 10],'electrodes','off');
    title([ 'Voltage, ' num2str(times2plot(i)) ' ms' ])
    
    % plot Laplacian map (spatially filtered)
    subplot(2,length(times2plot),i+length(times2plot))
    topoplot(laplacian_perrinX(tempdata,X,Y,Z),EEG.chanlocs,'plotrad',.53,'maplimits',[-40 40],'electrodes','off');
    title([ 'Lap., ' num2str(times2plot(i)) ' ms' ])
end

%% 
