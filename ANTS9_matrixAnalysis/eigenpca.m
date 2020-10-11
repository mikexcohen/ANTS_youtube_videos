%% mikexcohen.com

%% eigendecomposition

A = [3 1; 1 2];
[eigvecs,eigvals] = eig(A);

% define vectors
v1 = [.7 -.5]';
v2 = eigvecs(:,1);

% matrix-vector multiplication
v1A = A*v1;
v2A = A*v2;

% maximum value for plotting
xval = max([ abs(v1A); abs(v2A) ])*1.1;


figure(1), clf

subplot(131)
imagesc(A), axis square
title('Matrix A')


subplot(132)
plot([0 v1(1)],[0 v1(2)],'k','linew',4)
hold on
plot([0 v1A(1)],[0 v1A(2)],'r--','linew',2)
axis square, axis([-xval xval -xval xval])
plot(get(gca,'xlim'),[0 0],'k:')
plot([0 0],get(gca,'ylim'),'k:')
legend({'v';'Av'})
title('Not an eigenvector!')



subplot(133)
plot([0 v2(1)],[0 v2(2)],'k','linew',4)
hold on
plot([0 v2A(1)],[0 v2A(2)],'r--','linew',2)
axis square, axis([-xval xval -xval xval])
plot(get(gca,'xlim'),[0 0],'k:')
plot([0 0],get(gca,'ylim'),'k:')
legend({'v';'Av'})
title('Yes an eigenvector!')

%% PCA on simulated data

% data
x = [ 1*randn(1000,1) .4*randn(1000,1) ];

% rotation matrix
th = pi/4;
R1 = [ cos(th) -sin(th);
       sin(th)  cos(th) ];

% rotate data
y = x*R1;

% PCA of y (correlated data)
y = bsxfun(@minus,y,mean(y,1));
covmat = (y'*y) / length(y);
[evecsY,evalsY] = eig(covmat);


figure(2), clf
subplot(121)
plot(y(:,1),y(:,2),'m.','markersize',5)
set(gca,'xlim',[-5 5],'ylim',[-5 5])
hold on
plot(evalsY(1,1)*[0 evecsY(1,1)],evalsY(1,1)*[0 evecsY(2,1)],'k','linew',4)
plot(evalsY(2,2)*[0 evecsY(1,2)],evalsY(2,2)*[0 evecsY(2,2)],'k','linew',4)
xlabel('x-axis'), ylabel('y-axis')
axis square

% compute component scores
pc1 = y*evecsY(:,1);
pc2 = y*evecsY(:,2);

subplot(122)
plot(pc2,pc1,'m.')
set(gca,'xlim',[-5 5],'ylim',[-5 5])
xlabel('PC1 axis'), ylabel('PC2 axis')
axis square

%% vectors vs. values

x = rand(20);
covmat = (x'*x) / length(x);
[evecsX,evalsX] = eig(covmat);

figure(3), clf
subplot(121)
imagesc(evecsX)
axis square

subplot(122)
imagesc(evalsX)
set(gca,'clim',[-.2 .2])
axis square

%% PCA on EEG data

load sampleEEGdata

% compute ERP
erp = squeeze(mean(EEG.data,3));

% subtract mean and compute covariance
erp      = bsxfun(@minus,erp,mean(erp,2));
covarERP = (erp*erp')./(EEG.pnts-1);

figure(4), clf
subplot(121)
imagesc(covarERP)
axis square
set(gca,'xticklabel',{EEG.chanlocs(get(gca,'xtick')).labels},'yticklabel',{EEG.chanlocs(get(gca,'ytick')).labels},'clim',[-1 5])
title('Covariance of ERP')


% average single-trial covariances
subplot(122)
covarST = zeros(EEG.nbchan);
% note that the covariance of each trial is computed separately, then averaged
for i=1:EEG.trials
    eeg     = bsxfun(@minus,squeeze(EEG.data(:,:,i)),squeeze(mean(EEG.data(:,:,i),2)));
    covarST = covarST + (eeg*eeg')./(EEG.pnts-1);
end
covarST = covarST./i;
imagesc(covarST)
axis square
set(gca,'xticklabel',{EEG.chanlocs(get(gca,'xtick')).labels},'yticklabel',{EEG.chanlocs(get(gca,'ytick')).labels},'clim',[20 150])
title('Average covariance of single-trial EEG')

%%

% PCA of ERP
[evecsERP,evalsERP] = eig(covarERP);
[evecsST ,evalsST]  = eig(covarST);

figure(5), clf
subplot(221)
topoplot(evecsERP(:,end),EEG.chanlocs,'maplimits',[0 .2]);
title('PCA of ERP, 1^s^t component')

subplot(222)
topoplot(evecsST(:,end),EEG.chanlocs,'maplimits',[0 .2]);
title('PCA of EEG, 1^s^t component')

%% now get time courses

% for ERP
erpV1 = erp'*evecsERP(:,end);

% for EEG
eeg   = reshape(EEG.data,EEG.nbchan,[]);
eegV1 = eeg'*evecsST(:,end);
eegV1 = reshape(eegV1',EEG.pnts,EEG.trials);
eegV1 = mean(eegV1,2);

subplot(223)
plot(EEG.times,erpV1,'linew',2)
set(gca,'xlim',[-200 1200],'ylim',[-50 50])

subplot(224)
plot(EEG.times,eegV1,'linew',2)
set(gca,'xlim',[-200 1200],'ylim',[-50 50])

%%


%% ICA vs. PCA

% generate data

% data
x = [ 1*randn(1000,1) .05*randn(1000,1) ];

% rotation matrix
th = -pi/6;
R1 = [ cos(th) -sin(th); sin(th) cos(th) ];
th = -pi/3;
R2 = [ cos(th) -sin(th); sin(th) cos(th) ];

% rotate data
y = [ x*R1 ; x*R2 ];


figure(6), clf
subplot(221)
plot(y(:,1),y(:,2),'o')

datarange = max(y(:))*1.2;
set(gca,'xlim',[-datarange datarange],'ylim',[-datarange datarange])
xlabel('X axis'), ylabel('Y axis')
axis square
title('Data in XY space')


% now PCA
y = bsxfun(@minus,y,mean(y,1));
covmat = (y'*y) / length(y);
[evecsY,evalsY] = eig(covmat);

hold on
plot(evalsY(1,1)*[0 evecsY(1,1)],evalsY(1,1)*[0 evecsY(2,1)],'r','linew',4)
plot(evalsY(2,2)*[0 evecsY(1,2)],evalsY(2,2)*[0 evecsY(2,2)],'r','linew',4)


subplot(222)
pc1 = y*evecsY(:,1);
pc2 = y*evecsY(:,2);

plot(pc2,pc1,'ms')
datarange = max([pc1(:); pc2(:)])*1.2;
set(gca,'xlim',[-datarange datarange],'ylim',[-datarange datarange])
xlabel('PC1 axis'), ylabel('PC2 axis')
axis square
title('Data in PC space')





% now ICA
subplot(223)
plot(y(:,1),y(:,2),'o')
datarange = max(y(:))*1.2;
set(gca,'xlim',[-datarange datarange],'ylim',[-datarange datarange])

ivecs = jader(y');
hold on
plot([0 ivecs(1,1)],[0 ivecs(2,1)],'r','linew',4)
plot([0 ivecs(1,2)],[0 ivecs(2,2)],'r','linew',4)
xlabel('X axis'), ylabel('Y axis')
axis square
title('Data in XY space')



subplot(224)
ic_scores = ivecs*y';
plot(ic_scores(1,:),ic_scores(2,:),'ms')
datarange = max(ic_scores(:))*1.2;
set(gca,'xlim',[-datarange datarange],'ylim',[-datarange datarange])
xlabel('IC1 axis'), ylabel('IC2 axis')
axis square
title('Data in IC space')


%%


