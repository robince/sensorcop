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
d2dat = zeros(size(fdat));
for eli=1:2
    for trli=1:size(fdat,1)
        ddat(trli,:,eli) = gradient(squeeze(fdat(trli,:,eli)),2);
        d2dat(trli,:,eli) = gradient(ddat(trli,:,eli),2);
    end
end

% resample
tmp = permute(fdat,[2 1 3]);
rsfdat = resample(tmp,1,4);
rsfdat = reshape(rsfdat,[ size(rsfdat,1) size(tmp,2) size(tmp,3)]);
rsfdat = permute(rsfdat,[2 1 3]);

tmp = permute(ddat,[2 1 3]);
rsddat = resample(tmp,1,4);
rsddat = reshape(rsddat,[ size(rsddat,1) size(tmp,2) size(tmp,3)]);
rsddat = permute(rsddat,[2 1 3]);

time = resample(time,1,4);

% tmp = permute(d2dat,[2 1 3]);
% rsd2dat = resample(tmp,1,4);
% rsd2dat = reshape(rsd2dat,[ size(rsddat,1) size(tmp,2) size(tmp,3)]);
% rsd2dat = permute(rsd2dat,[2 1 3]);

cstm = copnorm(eyestm);

% cdat = copnorm(fdat);
% cddat = copnorm(ddat);

cdat = copnorm(rsfdat);
cddat = copnorm(rsddat);

[Ntrl, Nt, ~] = size(cdat);
stmi = 1;
otherstm = 2;
eleci = 2;

%% full p,vinfo
Ipv = zeros(Nt,1);
for ti=1:Nt
    Ipv(ti) = cmi_ggg(cstm(:,stmi), [cdat(:,ti,eleci) cddat(:,ti,eleci) ], cstm(:,otherstm), true, false, false);
end

%%
%% full p,v info - cond lag 1
Ipvacc = zeros(Nt,1);
for ti=2:Nt
    Ipvacc(ti) = cmi_ggg(cstm(:,stmi), [cdat(:,ti,eleci) cddat(:,ti,eleci)], [cdat(:,ti-1,eleci) cddat(:,ti-1,eleci)  cstm(:,otherstm)], true, false, false);
end


%%
figure
subplot(211)
plot(time,Ipv)
xlim(xl)
subplot(212)
plot(Ipvacc)
xlim(xl)

%% full p,v info - cond lag 1 and peaks
% peaks = [191 215];
peaks = [49 55];
Ipvaccpeak = zeros(Nt,1);
for ti=2:Nt
    conddat = [cdat(:,ti-1,eleci) cddat(:,ti-1,eleci)  cstm(:,otherstm)];
    for pi=1:length(peaks)
        if ti>(peaks(pi)+1)
            conddat = [conddat cdat(:,peaks(pi),eleci) cddat(:,peaks(pi),eleci)];
        end
    end
    Ipvaccpeak(ti) = cmi_ggg(cstm(:,stmi), [cdat(:,ti,eleci) cddat(:,ti,eleci)], conddat, true, false, false);
end

%%
figure
subplot(311)
plot(time,Ipv)
xlim(xl)
subplot(312)
plot(time,Ipvacc)
xlim(xl)
subplot(313)
plot(time,Ipvaccpeak)
xlim(xl)

%%
figure
ax = [];
ax(1) = subplot(211);
plot(time,Ipv,'k','LineWidth',2)
xlim(xl)

ax(2) = subplot(212);
plot(time,Ipvaccpeak, 'k', 'LineWidth',2)
xlim(xl)

set(ax,'xgrid','on');