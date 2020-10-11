%% mikexcohen.com

%%

activity = 1:.01:20; % activity
baseline = 10; % baseline

db = 10*log10(activity./ baseline );
pc = 100*( activity-baseline)./ baseline;


figure(1), clf
plot(activity,'linew',3)
hold on
plot(db,'r','linew',3)
legend({'"activity"','dB'})

%% comparing dB and percent change

figure(2), clf
plot(db,pc)
xlabel('dB'), ylabel('Percent change')

% find indices where db is closest to -/+2
[~,dbOfplus2]  = min(abs(db-+2));
[~,dbOfminus2] = min(abs(db--2));


hold on
axislim=axis;
plot([db(dbOfplus2)  db(dbOfplus2)], [pc(dbOfplus2)  axislim(3)],'k',[axislim(1) db(dbOfplus2)], [pc(dbOfplus2)  pc(dbOfplus2)], 'k')
plot([db(dbOfminus2) db(dbOfminus2)],[pc(dbOfminus2) axislim(3)],'k',[axislim(1) db(dbOfminus2)],[pc(dbOfminus2) pc(dbOfminus2)],'k')

%%