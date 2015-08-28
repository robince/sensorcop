% compare bias between different info methods
subid = 'HKI1';
fname = sprintf('%s_ridat.mat',subid);
dat = load(fullfile(data_dir, fname));

time = dat.time;
stim = dat.stim;
% rspdat = permute(dat.eegdat, [3 2 1]);
rspdat = permute(dat.csddat, [3 2 1]);
[Ntrl Nt Nchan] = size(rspdat);

tidx = (time>0) & (time<300);
time = time(tidx);
rspdat = rspdat(:,tidx,:);
[Ntrl Nt Nchan] = size(rspdat);

%%
figure
erp = squeeze(mean(rspdat,1));
erp0 = squeeze(mean(rspdat(stim==0,:,:),1));
erp1 = squeeze(mean(rspdat(stim==1,:,:),1));

chi = 41;
plot(time, erp0(:,chi),'r')
hold on
plot(time, erp1(:,chi),'k')
xlabel('Time (ms)')
legend('Face','Noise','Location','NorthWest')
set(gcf,'Pos',[ 680   957   256   141])
%%
Ntrl = 30;
chi = 41;
% chi = 24;
Nthread = 4;
Nt = length(time);

Ntrlstim = floor(Ntrl/2);
rsp0 = rspdat(stim==0,:,chi);
rsp1 = rspdat(stim==1,:,chi);

thsdat = [rsp0(randperm(size(rsp0,1),Ntrlstim),:); rsp1(randperm(size(rsp1,1),Ntrlstim),:)];
thsstm = [zeros(Ntrlstim,1); ones(Ntrlstim,1)];

cdat = copnorm(thsdat);
qstm = int16(thsstm);

Icop = info_c1d_slice_nobc_omp(cdat, qstm+1, 2, Ntrl, Nthread);

qdat = int16(bin.eqpop_slice(thsdat,2));
Ib2 = info.calc_info_slice_integer_c_int16_t(qdat,2,qstm,2,Ntrl);
qdat = int16(bin.eqpop_slice(thsdat,4));
Ib4 = info.calc_info_slice_integer_c_int16_t(qdat,4,qstm,2,Ntrl);
qdat = int16(bin.eqpop_slice(thsdat,8));
Ib8 = info.calc_info_slice_integer_c_int16_t(qdat,8,qstm,2,Ntrl);

t = -fastt2(thsdat, thsstm);
ks = kstest_slice_omp(thsdat, qstm, Ntrl, Nthread);

%%


figure
ax = [];
ax(1) = subplot(411);
plot(time,Icop);
axis tight
subplot(412)
plot(time,Ib2);
hold all
plot(time,Ib4);
plot(time,Ib8);
axis tight
set(ax(1),'YLim',get(gca,'YLim'));
subplot(413)
plot(time,t)
axis tight
subplot(414)
plot(time,ks)
axis tight
set(gcf,'Pos',[680   832   207   266])

%%
%% 2d WITH GRADIENT
Ntrl = 1000;
chi = 41;
% chi = 24;
Nthread = 4;
Nt = length(time);

Ntrlstim = floor(Ntrl/2);
rsp0 = rspdat(stim==0,:,chi);
rsp1 = rspdat(stim==1,:,chi);

rsp0 = cat(2, reshape(rsp0, size(rsp0,1), 1, []), reshape(gradient_dim1(rsp0')', size(rsp0,1), 1, []));
rsp1 = cat(2, reshape(rsp1, size(rsp1,1), 1, []), reshape(gradient_dim1(rsp1')', size(rsp1,1), 1, []));

thsdat = [rsp0(randperm(size(rsp0,1),Ntrlstim),:,:); rsp1(randperm(size(rsp1,1),Ntrlstim),:,:)];
thsstm = [zeros(Ntrlstim,1); ones(Ntrlstim,1)];

cdat = copnorm(thsdat);
qstm = int16(thsstm);

Icop = info_c1d_slice_nobc_omp(squeeze(cdat(:,1,:)), qstm+1, 2, Ntrl, Nthread);
Icop2 = info_cd_slice_nobc_omp(permute(cdat,[2 1 3]),2, qstm+1, 2, Ntrl, Nthread);

%%
figure
plot(time,Icop)
hold all
plot(time,Icop2)
axis tight
legend('V', '[V, dV/dt]')

%% temporal interaction
Ntrl = size(rspdat,1);
chi = 41;
% chi = 24;
Nthread = 4;
Nt = length(time);

thsdat = rspdat(:,:,chi);
thsstm = stim;

cdat = copnorm(thsdat);
qstm = int16(thsstm);

Icop = info_c1d_slice_nobc_omp(cdat, qstm+1, 2, Ntrl, Nthread);

%%
Iint = zeros(Nt,Nt);
for t1=1:Nt
    for t2=1:Nt
        I1 = info_gd(cdat(:,t1), stim, 2, false,false, false);
        if t1==t2
            Iint(t1,t2) = I1;
            continue
        end
        I2 = info_gd(cdat(:,t2), stim, 2, false,false, false);
        IJ = info_gd([cdat(:,t1) cdat(:,t2)],stim, 2,false, false, false);
        Iint(t1,t2) = I1 + I2 - IJ;
    end
end

% %%
% Iintn = zeros(Nt,Nt);
% for t1=1:Nt
%     for t2=1:Nt
%         Iintn(t1,t2) = 100*Iint(t1,t2) ./ max([Iint(t1,t1),Iint(t2,t2)]);
%     end
% end
% 
% imagesc(Iintn);colorbar

%%
rgb= [ ...
    94    79   162
    50   136   189
   102   194   165
   171   221   164
   230   245   152
   255   255   191
   254   224   139
   253   174    97
   244   109    67
   213    62    79
   158     1    66  ] / 255;

figure
imagesc(time,time,Iint);
colormap(intcm)
cx = caxis;
cxm = max(abs(cx));
caxis([-cxm cxm]);
cb = colorbar;
set(cb,'YLim', cx)
axis image
% box off

