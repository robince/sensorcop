subi = 1;
flti = 1;

subid = subjects{subi};
fname = sprintf('temporal_info_MIelec_%s_%s.mat',subid,fqnames{flti});
subres = matfile(fullfile(data_dir,fname));

fname = sprintf('data_%s_%s.mat',subid,fqnames{flti});
dat = matfile(fullfile(data_dir,fname));
eyestm = dat.eyebubs;
eegdat = dat.ferp(:,:,dat.RE);

    
Nbin = 10;
qeyestm = bin.eqpop_slice(eyestm,Nbin);
leye = 1;
reye = 2;
xl = [0 400];

t = subres.p;
Ip = t.Itc(:,1,2);
xIp = t.intItc(:,:,1,2);
t = subres.pv;
Ipv = t.Itc(:,1,2);
xIpv = t.intItc(:,:,1,2);
[Ntrl Nt] = size(eegdat);
rc = zeros(Nt,1);
for ti=1:Nt
    rc(ti) = corr(eyestm(:,1), eegdat(:,ti), 'type', 'spearman');
end
%%
figure
hold all
plot(time,Ip)
plot(time,Ipv)


%% calculate gradient info
subi = 1;
flti = 1;

subid = subjects{subi};
fname = sprintf('data_%s_%s.mat',subid,fqnames{flti});
dat = matfile(fullfile(data_dir,fname));
eyestm = dat.eyebubs;
otl = dat.ferp(:,:,dat.LEmaxMI);
otr = dat.ferp(:,:,dat.REmaxMI);
fdat = cat(3, otl, otr);

ddat = zeros(size(fdat));
for eli=1:2
    for trli=1:size(fdat,1)
        ddat(trli,:,eli) = gradient(squeeze(fdat(trli,:,eli)),2);
    end
end


cstm = copnorm(eyestm);
cddat = copnorm(ddat);

[Ntrl, Nt, ~] = size(cddat);
stmi = 1;
otherstm = 2;
eleci = 2;

Iv = zeros(Nt,1);
for ti=1:Nt
    Iv(ti) = cmi_ggg(cstm(:,stmi), cddat(:,ti,eleci), cstm(:,otherstm), true, false, false);
end

%%
figure
ax = [];
ax(1) = subplot(5,1,1);
% plot(time, erpdat.p.erp(:,el), 'k', 'LineWidth',2)
hold on
cidx = round(linspace(1, 32, Nbin));
for bi=0:(Nbin-1)
    thsbin = squeeze(mean(eegdat(qeyestm(:,leye)==(Nbin-bi-1), :),1));
    plot(time, thsbin, 'color', blue(cidx(end-bi),:), 'LineWidth', 2);
end
box on

ax(2) = subplot(5,1,2);
plot(time, (rc), 'k', 'LineWidth',2)
hline(0,':k')

yl = [-0.05 0.25];
ax(3) = subplot(5,1,3);
plot(time, Ip, 'k', 'LineWidth',2)
ylim(yl)
yl = [-0.05 0.25];

ax(4) = subplot(5,1,4);
plot(time, Iv, 'k', 'LineWidth',2)
ylim(yl)


ax(5) = subplot(5,1,5);
plot(time, Ip, 'LineWidth', 2, 'color', [0.8 0.8 0.8])
hold on
plot(time, Iv, 'LineWidth', 2, 'color', [0.8 0.8 0.8])
plot(time, Ipv, 'k', 'LineWidth', 2)
ylim(yl)

set(ax,'XLim',xl)
set(ax(1),'Color',[0.8 0.8 0.8]);
set(ax,'xgrid','on');

%%
figure
subplot(121)
imagesc(xIp)
cl = caxis;
cl = max(abs(cl));
cl = [-cl cl];
caxis(cl)
colorbar
axis square

subplot(122)
imagesc(xIpv)
cl = caxis;
cl = max(abs(cl));
cl = [-cl cl];
caxis(cl)
colorbar
axis square

%%
xl = [0 400];
idx = (time>xl(1)) & (time<xl(2));
intNt = sum(idx);
inttime = time(idx);

figure
axm = subplot(5,5,[2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20 ]);

% Iint = -res.(['int' Ivar])(:,:,2,1);
% I = res.(Ivar)(idx,2,1);
% II = res.intIelec(:,:,1);
% imagesc(inttime,inttime,-Iint);
% nint = normalise_interaction_matrix(-Iint,I,I,II,subthr);
imagesc(inttime,inttime,xIp);
% imagesc(nint);
% caxis([0 100])
cl = caxis;
cl = max(abs(cl));
cl = [-cl cl];
caxis(cl)
colorbar
axis square
% colormap pmkmp
colormap parula
title('Temporal interaction')
% caxis([-0.02 0])
cl = caxis;
% caxis([min(cl) 0])

lw = 2;
subplot(5,5,[22 23 24 25])
plot(time, Ip, 'k', 'LineWidth', lw);
axis tight
% set(gca,'Color',0.8*[1 1 1])
xlim(xl)
% tax = colorbar;
% set(tax,'Vis','off')
intpos = get(axm,'Pos');
pos = get(gca,'Pos');
pos(1) = intpos(1);
pos(3) = intpos(3);
set(gca,'Pos',pos)


box off

subplot(5,5,[1 6 11 16 ])
plot(time, Ip, 'k', 'LineWidth', lw);
axis tight
box off
xlim(xl)
% set(gca,'Color',0.8*[1 1 1])
set(gca,'CameraUpVector', [-1 0 0])
% set(gcf,'Position',[46         175        1012         782]);
% suptitle('Left eye')

%%
xl = [0 400];
idx = (time>xl(1)) & (time<xl(2));
intNt = sum(idx);
inttime = time(idx);

figure
axm = subplot(5,5,[2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20 ]);

% Iint = -res.(['int' Ivar])(:,:,2,1);
% I = res.(Ivar)(idx,2,1);
% II = res.intIelec(:,:,1);
% imagesc(inttime,inttime,-Iint);
% nint = normalise_interaction_matrix(-Iint,I,I,II,subthr);
imagesc(inttime,inttime,-xIpv);
% imagesc(nint);
% caxis([0 100])
cl = caxis;
cl = max(abs(cl));
cl = [0 cl];
caxis(cl)
colorbar
axis square
% colormap pmkmp
colormap parula
title('Temporal interaction')
% caxis([-0.02 0])
cl = caxis;
% caxis([min(cl) 0])

lw = 2;
subplot(5,5,[22 23 24 25])
plot(time, Ipv, 'k', 'LineWidth', lw);
axis tight
% set(gca,'Color',0.8*[1 1 1])
xlim(xl)
% tax = colorbar;
% set(tax,'Vis','off')
intpos = get(axm,'Pos');
pos = get(gca,'Pos');
pos(1) = intpos(1);
pos(3) = intpos(3);
set(gca,'Pos',pos)


box off

subplot(5,5,[1 6 11 16 ])
plot(time, Ipv, 'k', 'LineWidth', lw);
axis tight
box off
xlim(xl)
% set(gca,'Color',0.8*[1 1 1])
set(gca,'CameraUpVector', [-1 0 0])
% set(gcf,'Position',[46         175        1012         782]);
% suptitle('Left eye')