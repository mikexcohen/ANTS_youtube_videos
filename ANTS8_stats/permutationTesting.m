% mikexcohen@gmail.com

% load in data
load v1_laminar.mat

% useful variables for later...
npnts   = size(csd,2);
ntrials = size(csd,3);

%% setup some stuff and junk

% pick channels
% 1 is deep (hippocampus), ~7 is Layer-4 in V1
% hint: try with other channel pairs!
chan1idx = 9;
chan2idx = 5;

times2save = -.3:.01:1.1; % time vector is in seconds!
xlim       = [-.1 1]; % for plotting

% specify frequency range
min_freq = 10;
max_freq = 100;
num_frex = 50;

% define frequency and time ranges
frex = linspace(
times2saveidx = dsearchn(timevec',times2save');

% parameters for complex Morlet wavelets
wavtime  = -2:1/srate:2-1/srate;
half_wav = (length(wavtime)-1)/2;
cycRange = [ 4 10 ];
nCycles  = logspace(log10(cycRange(1)),log10(cycRange(end)),num_frex);

% FFT parameters
nWave = length(wavtime);
nData = npnts*ntrials;
nConv = nWave+nData-1;

% and create wavelets
cmwX = zeros(num_frex,nConv);
for fi=1:num_frex
    s       =  % frequency-normalized width of Gaussian
    cmw      = exp(1i*2*pi*frex(fi).*wavtime) .* exp( (-wavtime.^2) ./ (2*s^2) );
    tempX     = fft(cmw,nConv);
    cmwX(fi,:) = tempX ./ max(tempX);
end

%% run convolution to extract tf power

% spectra of data
dataX1 = fft( reshape(csd(chan1idx,:,:),1,[]) ,nConv);
dataX2 = fft( reshape(csd(chan2idx,:,:),1,[]) ,nConv);


% initialize output time-frequency data
% notice that here we save all trials
tf = zeros(2,num_frex,length(times2save),ntrials);

for fi=1:num_frex
    
    % run convolution
    as1 = ifft(cmwX(fi,:).*dataX1);
    as1 = as1(half_wav+1:end-half_wav);
    as1 = reshape(as1,npnts,ntrials);
    
    % power on all trials from channel "1"
    % only from times2saveidx!
    tf(1,fi,:,:) = 
    
    
    % run convolution
    as2 = ifft(cmwX(fi,:).*dataX2);
    as2 = as2(half_wav+1:end-half_wav);
    as2 = reshape(as2,npnts,ntrials);
    
    % power on all trials from channel "2"
    % only from times2saveidx!
    tf(2,fi,:,:) = 
end

% for convenience, compute the difference in power between the two channels
diffmap = squeeze(mean(tf(2,:,:,:),4 )) - squeeze(mean(tf(1,:,:,:),4 ));


%% plotting the raw data

clim = [0 20000];

figure(1), clf
subplot(221)
imagesc(times2save,frex,squeeze(mean( tf(1,:,:,:),4 )))
set(gca,'clim',clim,'ydir','n','xlim',xlim)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Channel ' num2str(chan1idx) ])

subplot(222)
imagesc(times2save,frex,squeeze(mean( tf(2,:,:,:),4 )))
set(gca,'clim',clim,'ydir','n','xlim',xlim)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Channel ' num2str(chan2idx) ])

subplot(223)
imagesc(times2save,frex,diffmap)
set(gca,'clim',[-mean(clim) mean(clim)],'ydir','n','xlim',xlim)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Difference: channels ' num2str(chan2idx) ' - ' num2str(chan1idx) ])

%% statistics via permutation testing

% p-value
pval = 0.05;

% convert p-value to Z value
zval = abs(norminv(pval));

% number of permutations
n_permutes = 1000;

% initialize null hypothesis maps
permmaps = zeros(n_permutes,num_frex,length(times2save));

% for convenience, tf power maps are concatenated
%   in this matrix, trials 1:ntrials are from channel "1" 
%   and trials ntrials+1:end are from channel "2"
tf3d = cat(3,squeeze(tf(1,:,:,:)),squeeze(tf(2,:,:,:)));


% generate maps under the null hypothesis
for permi = 1:n_permutes
    
    % randomize trials, which also randomly assigns trials to channels
    randorder = randperm(size(tf3d,3));
    temp_tf3d = tf3d(:,:,randorder);
    
    % compute the "difference" map
    % what is the difference under the null hypothesis?
    permmaps(permi,:,:) = squeeze( mean(temp_tf3d(:,:,1:ntrials),3) - mean(temp_tf3d(:,:,ntrials+1:end),3) );
end

%% show non-corrected thresholded maps

% compute mean and standard deviation maps
mean_h0 = squeeze(mean(permmaps));
std_h0  = squeeze(std(permmaps));

% now threshold real data...
% first Z-score
zmap = (diffmap-mean_h0) ./ std_h0;

% threshold image at p-value, by setting subthreshold values to 0
zmap(abs(zmap)<zval) = 0;


% now some plotting...

figure(2), clf

subplot(221)
imagesc(times2save,frex,diffmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','nor')
title('TF map of real power values')

subplot(222)
imagesc(times2save,frex,diffmap);
hold on
contour(times2save,frex,logical(zmap),1,'linecolor','k');
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')
title('Power values and outlined significance regions')

subplot(223)
imagesc(times2save,frex,zmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-10 10],'xlim',xlim,'ydir','no')
title('Thresholded TF map of Z-values')

%% 






%% corrections for multiple comparisons

% initialize matrices for cluster-based correction
max_cluster_sizes = zeros(1,n_permutes);
% ... and for maximum-pixel based correction
max_val = zeros(n_permutes,2); % "2" for min/max

% loop through permutations
for permi = 1:n_permutes
    
    % take each permutation map, and transform to Z
    threshimg = squeeze(permmaps(permi,:,:));
    threshimg = (threshimg-mean_h0)./std_h0;
    
    % threshold image at p-value
    threshimg(abs(threshimg)<zval) = 0;
    
    
    % find clusters (need image processing toolbox for this!)
    islands = bwconncomp(threshimg);
    if numel(islands.PixelIdxList)>0
        
        % count sizes of clusters
        tempclustsizes = cellfun(@length,islands.PixelIdxList);
        
        % store size of biggest cluster
        max_cluster_sizes(permi) = max(tempclustsizes);
    end
    
    
    % get extreme values (smallest and largest)
    temp = sort( reshape(permmaps(permi,:,:),1,[] ));
    max_val(permi,:) = [ min(temp) max(temp) ];
    
end

%% show histograph of maximum cluster sizes

figure(3), clf
hist(max_cluster_sizes,20);
xlabel('Maximum cluster sizes'), ylabel('Number of observations')
title('Expected cluster sizes under the null hypothesis')


% find cluster threshold (need image processing toolbox for this!)
% based on p-value and null hypothesis distribution
cluster_thresh = prctile(max_cluster_sizes,100-(100*pval));

%% plots with multiple comparisons corrections

% now find clusters in the real thresholded zmap
% if they are "too small" set them to zero
islands = bwconncomp(zmap);
for i=1:islands.NumObjects
    % if real clusters are too small, remove them by setting to zero!
    if numel(islands.PixelIdxList{i}==i)<cluster_thresh
        zmap(islands.PixelIdxList{i})=0;
    end
end

% plot tresholded results
figure(4), clf
subplot(221)
imagesc(times2save,frex,diffmap)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('TF power, no thresholding') 
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')

subplot(222)
imagesc(times2save,frex,diffmap)
hold on
contour(times2save,frex,logical(zmap),1,'linecolor','k')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('TF power with contour')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')

subplot(223)
imagesc(times2save,frex,zmap)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('z-map, thresholded')
set(gca,'clim',[-13 13],'xlim',xlim,'ydir','normal')

%% now with max-pixel-based thresholding

% find the threshold for lower and upper values
thresh_lo = prctile(max_val(:,1),100*pval); % what is the
thresh_hi = prctile(max_val(:,2),100-100*pval); % true p-value?

% threshold real data
zmap = diffmap;
zmap(zmap>thresh_lo & zmap<thresh_hi) = 0;

figure(5), clf
subplot(221)
imagesc(times2save,frex,diffmap)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('tf power map, no thresholding') 
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','n')

subplot(222)
imagesc(times2save,frex,diffmap)
hold on
contour(times2save,frex,logical(zmap),1,'linecolor','k')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('tf power map with contour')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','normal')

subplot(223)
imagesc(times2save,frex,zmap)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('tf power map, thresholded')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','no')

%% end.
