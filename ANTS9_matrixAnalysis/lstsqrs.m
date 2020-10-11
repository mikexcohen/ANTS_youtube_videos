%% mikexcohen.com

%% matrix inverse

A = [2 3; 1 5];
Ainv = inv(A);

figure(1), clf
subplot(131)
imagesc(A)
axis off, axis square
title('A')

subplot(132)
imagesc(Ainv)
axis off, axis square
title('A^-^1')

subplot(133)
imagesc(A*Ainv)
axis off, axis square
title('AA^-^1 = I')

%% psuedo-inverse

A = [2 3; 1 3];
A = [2 3; 2 3];

inv(A)
pinv(A)

%%


%% least-squares: a simple example

n = 10; % number of data points

% generate some data
b = linspace(1,3,10)' + rand(n,1);

% design matrix for an intercept and linear effect
A = [ ones(n,1) (1:n)' ];

% compute the parameters
x = (A'*A)\A'*b;

% compute the model-predicted data 
yHat = x(1)*A(:,1) + x(2)*A(:,2);

% and plot!
figure(2), clf
plot(1:n,b,'o','markerface','b','markersize',20)
hold on
plot(1:n,yHat,'rp-','linew',2,'markersize',30,'markerface','k')
title('Correct')

set(gca,'xlim',[0 n+1],'ylim',[1 4.5])
legend({'observed data';'predicted data'})

%% load EEG data and extract reaction times in ms

load sampleEEGdata.mat

rts = zeros(size(EEG.epoch));

% loop over trials
for ei=1:EEG.trials
    
    % find the index corresponding to time=0, i.e., trial onset
    [~,zeroloc] = min(abs( cell2mat(EEG.epoch(ei).eventlatency) ));
    
    % reaction time is the event after the trial onset
    rts(ei) = EEG.epoch(ei).eventlatency{zeroloc+1};
end

% create design matrix
A = [ ones(EEG.trials,1) rts' ];

%% define convolution parameters for time-frequency analysis

freqrange = [2 20]; % extract only these frequencies (in Hz)
numfrex   = 30;     % number of frequencies between lowest and highest


% set up convolution parameters
wavtime = -2:1/EEG.srate:2;
frex    = linspace(freqrange(1),freqrange(2),numfrex);
nData   = EEG.pnts*EEG.trials;
nKern   = length(wavtime);
nConv   = nData + nKern - 1;
halfwav = (length(wavtime)-1)/2;
nCyc    = logspace(log10(4),log10(12),numfrex);

% create wavelets
cmwX = zeros(numfrex,nConv);
for fi=1:numfrex
    
    % create time-domain wavelet
    s   = nCyc(fi) / (2*pi*frex(fi));
    cmw = exp(2*1i*pi*frex(fi).*wavtime) .* exp( (-wavtime.^2) / (2*s.^2) );
    
    % compute fourier coefficients of wavelet and normalize
    cmwX(fi,:) = fft(cmw,nConv);
    cmwX(fi,:) = cmwX(fi,:) ./ max(cmwX(fi,:));
end


% initialize time-frequency output matrix
tf = zeros(numfrex,EEG.pnts);
tf3d = zeros(numfrex,EEG.pnts,EEG.trials);

% compute Fourier coefficients of EEG data (doesn't change over frequency!)
eegX = fft( reshape(EEG.data(47,:,:),1,[]) ,nConv);

% loop over frequencies
for fi=1:numfrex
    
    % second and third steps of convolution
    as = ifft( cmwX(fi,:).*eegX ,nConv );
    
    % cut wavelet back to size of data
    as = as(halfwav+1:end-halfwav);
    as = reshape(as,EEG.pnts,EEG.trials);
    
    % extract power from all trials
    tf3d(fi,:,:) = abs(as).^2;
    
end % end frequency loop

%% now compute correlations

% reshape the 3D matrix to 2D
tf2d = reshape(tf3d,numfrex*EEG.pnts,EEG.trials)';

% the 2D matrix can be used in a single least squares equation
x = (A'*A)\A'*tf2d;
covmat = reshape(x(2,:),numfrex,EEG.pnts); % demushing

%% show the design and data matrices

figure(3), clf

ax1_h = axes;
set(ax1_h,'Position',[.05 .1 .1 .8])
imagesc(A)
set(ax1_h,'xtick',1:2,'xticklabel',{'Int';'RTs'},'ydir','norm')
ylabel('Trials')


ax2_h = axes;
set(ax2_h,'Position',[.25 .1 .7 .8])
imagesc(tf2d)
set(ax2_h,'ydir','norm','clim',[0 20])
ylabel('Trials')
xlabel('Timefrequency')

colormap gray

%% show the results

figure(4), clf

% show time-frequency map of regressors
contourf(EEG.times,frex,covmat,40,'linecolor','none')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'xlim',[-200 1200],'clim',[-.012 .012])

%%
