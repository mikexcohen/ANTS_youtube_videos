%% group-level analyses
% mikexcohen.com

clear, close all

% Get a list of data files ready to be analyzed
sublist = dir('*tf.mat');
sublist = {sublist.name};

%% load in level-1 data

for subno=1:length(sublist)
    
    load(sublist{subno})
    
    % initialize matrices on 1st subject
    if subno==1
        tf_all = zeros([ length(sublist) size(mw_tf) ]);
        % list more variables here as applicable...
    end
    
    % place single-subject data into group-size matrices
    tf_all(subno,:,:,:,:,:) = mw_tf;
    
end % end loop around subjects

%% moving right along...

% This organization of the rest of the script is to make sections for each
% group-level analysis you want to do. At the start of each cell you can
% specify different parameters (e.g., different conditions or channels or
% time/frequency windows, etc.)

%% Plot time-frequency maps for specific channels

% plotting parameters
channel2plot = 'c4';

phsclim = [0 .4];
powclim = [-2 2];


chanidx = find(strcmpi(channel2plot,{chanlocs.labels}));

% same data, different plotting styles
figure(1), clf
for condi=1:3
    subplot(2,3,condi)
    contourf(tx,mw_frex,squeeze(mean(tf_all(:,condi,chanidx,:,:,1),1)),50,'linecolor','none')
    set(gca,'clim',powclim,'xlim',[-300 800],'yscale','log','ytick',ceil(logspace(0,log10(length(mw_frex)),8)));
    xlabel('Peri-stim time (ms)')
    ylabel('Frequency (Hz)')
    title([ chanlocs(chanidx).labels ' power; ' condition_labels{condi} ])
    
    subplot(2,3,condi+3);
	contourf(tx,mw_frex,squeeze(mean(tf_all(:,condi,chanidx,:,:,2),1)),50,'linecolor','none')
	set(gca,'clim',phsclim,'xlim',[-300 800],'yscale','log','ytick',ceil(logspace(0,log10(length(mw_frex)),8)));
    xlabel('Peri-stim time (ms)')
    ylabel('Frequency (Hz)')
    title([ chanlocs(chanidx).labels ' ITPC; ' condition_labels{condi} ])
end

%% topographical maps at specific time-frequency window

times2plot = [ 100 300 450 500];
freqrange  = [  8 12 ];
clim       = [ -2  2 ];

% convert times and frequencies to indices
timesidx = dsearchn(tx',times2plot');
freqsidx = dsearchn(mw_frex',freqrange');


figure(2), clf
set(gcf,'name','some nice topographical maps')
% loop over time points
for ti=1:length(times2plot)
    % loop over conditions
    for condi=1:3
        
        % specify place to plot
        subplot(3,length(times2plot),ti+(condi-1)*length(times2plot))
        
        % make topomap
		topoplot(squeeze(mean(mean(tf_all(:,condi,:,freqsidx(1):freqsidx(2),timesidx(ti),1),1),4)),chanlocs,'plotrad',.55,'maplimits',clim,'electrodes','off','numcontour',0);
        title([ condition_labels{condi} ' (' num2str(times2plot(ti)) 'ms; ' num2str(freqrange(1)) '-' num2str(freqrange(2)) 'Hz)' ]);
    end
end

%% Extract statistics: One window for all subjects/conditions
% Here, we will extract data from time-frequency windows of interest 
% and put them into a text file that can easily be read by SPSS.

time_windows = [ 0 200; 200 400 ]; % two time windows, boundaries separated by semicolon
freq_windows = [ 4 8 ]; % of course you could have multiple frequency windows
chans4analysis={'C3';'c4';'po7';'po8'};

statsfilename='statistics_file.txt';


 %%%%%%%%%

% find indices corresponding to time and frequency windows, and electrodes
timeidx  = zeros(size(time_windows));
freqidx  = zeros(size(freq_windows));
chansidx = zeros(size(chans4analysis));

for i=1:size(time_windows,1)
    for j=1:2
        [~,timeidx(i,j)] = min(abs(tx-time_windows(i,j)));
    end
end
for i=1:size(freq_windows,1)
    for j=1:2
        [~,freqidx(i,j)] = min(abs(mw_frex-freq_windows(i,j)));
    end
end
for i=1:length(chans4analysis)
    chansidx(i) = find(strcmpi(chans4analysis{i},{chanlocs.labels}));
end


% pointer to stats file
fid=fopen(statsfilename,'w');

% now write out data
for subno=0:size(tf_all,1)
    
    % Write out column variable name or subject number.
    if subno==0
        fprintf(fid,'subnum\t');
    else
        fprintf(fid,'%g\t',subno);
    end
    
    for chani=1:length(chans4analysis)
        for condi=1:length(condition_labels)
            for ti=1:size(time_windows,1)
                for fi=1:size(freq_windows,1)
                    
                    % Write out the column variable name of data.
                    if subno==0
                        fprintf(fid,[ chanlocs(chansidx(chani)).labels '_' condition_labels{condi} '_t' num2str(time_windows(ti,1)) num2str(time_windows(ti,2)) '_f' num2str(freq_windows(fi,1)) num2str(freq_windows(fi,2)) '\t' ]);
                    else
                        fprintf(fid,'%g\t',mean(mean(tf_all(subno,condi,chansidx(chani),freqidx(fi,1):freqidx(fi,2),timeidx(ti,1):timeidx(ti,2),1),4),5));
                    end
                    
                end % end frequency loop
            end % end time loop
        end % end condition loop
    end % end channel loop
    fprintf(fid,'\n');
end % end subject loop

fclose(fid);

%% Constrained subject-specific TF windows.

time_windows = [ 200 800 ];
freq_windows = [ 2 10 ];
chans4analysis={'C3';'c4';'po7';'po8'};

statsfilename='statistics_file_subjectSpecific.txt';


 %%%%%%%%%
 
% find indices corresponding to time and frequency windows, and electrodes
timeidx=zeros(size(time_windows));
freqidx=zeros(size(freq_windows));
chansidx=zeros(size(chans4analysis));

for i=1:length(time_windows)
    [junk,timeidx(i)]=min(abs(tx-time_windows(i)));
end
for i=1:length(freq_windows)
    [junk,freqidx(i)]=min(abs(mw_frex-freq_windows(i)));
end
for i=1:length(chans4analysis)
    chansidx(i)=find(strcmpi(chans4analysis{i},{chanlocs.labels}));
end


% pointer to stats file
fid=fopen(statsfilename,'w');

% now write out data
for subno=0:size(tf_all,1)
    
    % Write out column variable name or subject number.
    if subno==0
        fprintf(fid,'subnum\t');
    else
        fprintf(fid,'%g\t',subno);
    end
    
    for chani=1:length(chans4analysis)
        for condi=1:length(condition_labels)
            
            % Write out the column variable name of data.
            if subno==0
                fprintf(fid,[ chans4analysis{chani} '_' condition_labels{condi} '_t' num2str(time_windows(1)) num2str(time_windows(2)) '_f' num2str(freq_windows(1)) num2str(freq_windows(2)) '\t' ]);
                fprintf(fid,[ chans4analysis{chani} '_' condition_labels{condi} '_t' num2str(time_windows(1)) num2str(time_windows(2)) '_f' num2str(freq_windows(1)) num2str(freq_windows(2)) '_freq\t' ]);
                fprintf(fid,[ chans4analysis{chani} '_' condition_labels{condi} '_t' num2str(time_windows(1)) num2str(time_windows(2)) '_f' num2str(freq_windows(1)) num2str(freq_windows(2)) '_time\t' ]);
            else
                
                % get data from large window
                % Note: This is condition-specific; you could make it condition-average by taking the mean across the 2nd dimension.
                tempdat = squeeze(tf_all(subno,condi,chansidx(chani),freqidx(1):freqidx(2),timeidx(1):timeidx(2),1));
                
                % find max point (must first vectorize otherwise you get max for each line)
                [junk,maxpoint]=max(tempdat(:));
                
                % convert that back to 2d coordinates
                [maxF,maxT]=ind2sub(size(tempdat),maxpoint);
                
                % get the indices from the larger matrix, not the indices of the smaller TF window.
                maxF = maxF+freqidx(1)-1;
                maxT = maxT+timeidx(1)-1;
                
                % Now write out data and max time/frequency to file.
                fprintf(fid,'%g\t',mean(mean(tf_all(subno,condi,chansidx(chani),maxF-1:maxF+1,maxT-1:maxT+1,1),4),5));
                fprintf(fid,'%g\t',mw_frex(maxF));
                fprintf(fid,'%g\t',tx(maxT));
                
                % this could also be one line:
                % fprintf(fid,'%g\t%g\t%g\t',mean(mean(tf_all(subno,condi,chani,maxF-1:maxF+1,maxT-1:maxT+1,1),4),5),mw_frex(maxF),tx(maxT));
            end
        end % end condition loop
    end % end channel loop
    fprintf(fid,'\n');
end % end subject loop

fclose(fid);

%% end.

