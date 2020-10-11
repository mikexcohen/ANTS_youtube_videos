% eeglab's EEG structure and indexing
% mikeXcohen@gmail.com

%% load in EEG data

load sampleEEGdata.mat

% FYI, this would also work:
% [file2load,path4file]=uigetfile('*.mat','Please select EEG data file');
% load([ path4file file2load ])

% take a minute to inspect the EEG structure
EEG

%% finding time indices based on ms

% The problem: We want to create a topographical plot at time=300 ms

time2plot = 300; % in ms!

% extract the trial-averaged data from requested time point
% (note: this code contains an error!)
data2plot = squeeze(mean( EEG.data(:,time2plot,:) ,3));


% plot
figure(1), clf
topoplot(data2plot,EEG.chanlocs);
title([ 'Topoplot from time=' num2str(EEG.times(time2plot)) ' ms.' ])

%% same concept for frequencies

frex = linspace(2,100,42);

freqIwant = 23; % in hz

% use min(abs trick to find closest frequency to 23 Hz
[~,frexidx] = min(abs(frex-freqIwant));

% the function dsearchn also works
frexidx = dsearchn(frex',freqIwant);

%% indexing channels based on names

% the electrode label that we want to analyze
electrodeName = 'p1'; % case doesn't matter

% find the channel number that corresponds to this label
electrodeidx = strcmpi(electrodeName,{EEG.chanlocs.labels});

% confirm that this electrode is correct
EEG.chanlocs(electrodeidx)

% plot the ERP from this electrode
figure(1), clf
plot(EEG.times,mean( EEG.data(electrodeidx,:,:),3 ))


%% now multiple electrodes

electrodeNames = {'p1';'fc6';'t18'};

% initialize
electrodeidx = zeros(size(electrodeNames));

% loop through electrodes and find the index of each one
% (The code below has three errors. Try to find/fix them!)
for chani=1:length(electrodeNames)
    electrodeidx(chani) = strcmpi(electrodeNames,{EEG.chanlocs.labels});
end

% plot all the ERPs
plot(EEG.times,mean( EEG.data(electrodeidx,:,:),3 ))
legend({EEG.chanlocs(electrodeidx).labels})

%%

