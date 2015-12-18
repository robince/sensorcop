% load from sensor_info.m

chi = 40;
ti = 320;
% eeg = squeeze(dat.eegdat(chi,:,:));
csd = squeeze(dat.csddat(chi,:,:));
eeg = csd;

%%
Nplt = 10;
cols = {'r' 'b'};

figure
hold on
offset = 4;
pltidx = randperm(size(eeg,2),Nplt);
for pii=1:Nplt
    pi = pltidx(pii);
    plot(time, (pii-1)*offset + eeg(:,pi), cols{stim(pi)+1});
    plot(time(ti), (pii-1)*offset + eeg(ti,pi), '.','MarkerSize',16,'Color', cols{stim(pi)+1});
end
xlim([0 300])
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02 0.025])

%%
ens0 = eeg(ti,stim==0);
ens1 = eeg(ti,stim==1);

kxi = linspace(-3,3,100);

[k0,~] = ksdensity(ens0, kxi);
[k1,~] = ksdensity(ens1, kxi);

figure
plot(kxi,k0, cols{1},'LineWidth',2);
hold on
plot(kxi,k1, cols{2},'LineWidth',2);
box off
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02 0.025])