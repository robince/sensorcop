% calculate sensor info to pick a good channel

subid = 'HKI1';
fname = sprintf('%s_ridat.mat',subid);
dat = load(fullfile(data_dir, fname));

%%
time = dat.time;
stim = dat.stim;
rspdat = permute(dat.eegdat, [3 2 1]);
% rspdat = permute(dat.csddat, [3 2 1]);
[Ntrl Nt Nchan] = size(rspdat);

qdat = copnorm(rspdat);
%%
I = zeros(Nt, Nchan);
for ci=1:Nchan
    for ti=1:Nt
        I(ti,ci) = info_gd(qdat(:,ti,ci), stim, 2, true, true, false);
    end
end

%%
minlim = min(I(:));
maxlim = max(I(:));

tidx = find((time>0) & (time<350));

%%
clear F
for tii=1:length(tidx)
    ti = tidx(tii);
    t = I(ti,:);
    cla
    topoplot(t, dat.chanlocs, 'maplimits', [minlim maxlim]);
    title(sprintf('%d ms', time(ti)),'FontSize',12);
    F(tii) = getframe(gcf);
end

% topoplot(mean(DATA(:,170,:),3),EEG.locs)

fname = sprintf('csd_face_info.mov');
mo = QTWriter(fname);
mo.FrameRate = 5;
mo.Loop = 'loop';

for i=1:length(F)
    writeMovie(mo, F(i));
end
close(mo);
close all

%%
ci = 41; % chan
tidx = find((time>-100) & (time<350));
figure
subplot(311)
hold on
erp0 = mean(rspdat(stim==0,tidx,ci));
erp1 = mean(rspdat(stim==1,tidx,ci));
plot(time(tidx), erp0,'r')
plot(time(tidx), erp1,'b')
legend('Face','Noise')
title('ERPs')

subplot(312)
plot(time(tidx),erp1-erp0,'k');
title('ERP DIFF')

subplot(313)
plot(time(tidx),I(tidx,ci),'k');
title('Info')
