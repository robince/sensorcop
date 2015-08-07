subid = 'HKI1';
fname = sprintf('%s_ridat.mat',subid);
dat = load(fullfile(data_dir, fname));

time = dat.time;
stim = dat.stim;
% rspdat = permute(dat.eegdat, [3 2 1]);
rspdat = permute(dat.csddat, [3 2 1]);
[Ntrl Nt Nchan] = size(rspdat);

%%
tidx = (time>50) & (time<500);
time = time(tidx);
rspdat = rspdat(:,tidx,:);
[Ntrl Nt Nchan] = size(rspdat);

%%
Nthread = 16;

Nbin = 2;
qrsp = int16(bin.eqpop_slice_omp(rspdat(:,:),Nbin,Nthread));
Ib2 = reshape(info.calc_info_slice_omp_integer_c_int16_t(qrsp, Nbin, int16(stim), 2, Ntrl, Nthread), [Nt Nchan]);

Nperm = 1000;
Ib2perm = zeros(Nt, Nchan, Nperm);
for pi=1:Nperm
    pstim = int16(stim(randperm(Ntrl)));
    thsI = reshape(info.calc_info_slice_omp_integer_c_int16_t(qrsp, Nbin, pstim, 2, Ntrl, Nthread), [Nt Nchan]);
    Ib2perm(:,:,pi) = thsI;
end

maxMI = squeeze(max(max(Ib2perm,[],1),[],2));
thrMI = prctile(maxMI, 99);

%%
figure
subplot(121)
imagesc(time,[],Ib2')
colorbar
colormap parula
subplot(122)
imagesc(time,[],Ib2' > thrMI)
% colormap gray

%%
Ntrl = 500;
rspdat = rspdat(1:Ntrl,:,:);
stim = stim(1:Ntrl);
Nthread = 32;
Nbin = 16;
qrsp = copnorm_slice_omp_c_double(rspdat(:,:),Nthread);

% Icop = zeros([Nt Nchan]);
% for ci=1:size(qrsp,2);
%     Icop(ci) = info_gd( qrsp(:,ci), stim, 2, false, true, false );
% end
Icop = reshape(info_c1d_slice_nobc_omp(qrsp, int16(stim+1), 2, Ntrl, Nthread),[Nt Nchan]);

Nperm = 1000;
tic
Icopperm = zeros(Nt, Nchan, Nperm);
for pi=1:Nperm
    pstim = stim(randperm(Ntrl))+1;
    thsI = reshape(info_c1d_slice_nobc_omp(qrsp, int16(pstim), 2, Ntrl, Nthread),[Nt Nchan]);
    Icopperm(:,:,pi) = thsI;
end
toc

maxMI = squeeze(max(max(Icopperm,[],1),[],2));
thrMI = prctile(maxMI, 99);

%%
figure
subplot(121)
imagesc(time,[],Icop')
colorbar
colormap parula
subplot(122)
imagesc(time,[],Icop' > thrMI)
% colormap gray


%%
tic
t = reshape(fastt2(rspdat, stim),[Nt Nchan]);
toc
Nperm = 1000;
tic
tperm = zeros(Nt, Nchan, Nperm);
parfor pi=1:Nperm
    pstim = stim(randperm(Ntrl));
    thst = reshape(fastt2(rspdat, pstim),[Nt Nchan]);
    tperm(:,:,pi) = thst;
end
toc

maxt = squeeze(max(max(tperm,[],1),[],2));
mint = squeeze(min(min(tperm,[],1),[],2));
thrt_hi = prctile(maxt, 99.5);
thrt_lo = prctile(mint, 0.5);

%%

figure
subplot(121)
imagesc(time,[],t')
colorbar
colormap parula
subplot(122)
imagesc(time,[],(t' > thrt_hi) | (t'<thrt_lo))
% colormap gray

%%
tic
rsp = rspdat(:,:);
ks = zeros([Nt Nchan]);
parfor ci=1:size(rsp,2)
    stat = kstest2ri(rsp(stim==0,ci),rsp(stim==1,ci)); 
    ks(ci) = stat;
end
toc

Nperm = 1000;
tic
ksperm = zeros(Nt, Nchan, Nperm);
parfor pi=1:Nperm
    pstim = stim(randperm(Ntrl));
    thsks = zeros([Nt Nchan]);
    thsrsp1 = rsp(pstim==0,:);
    thsrsp2 = rsp(pstim==1,:);
    for ci=1:size(rsp,2)
        thsks(ci) = kstest2ri(thsrsp1(:,ci),thsrsp2(:,ci));
    end
    ksperm(:,:,pi) = thsks;
end
toc

maxks = squeeze(max(max(ksperm,[],1),[],2));
thrks = prctile(maxks, 99);

%%

figure
subplot(121)
imagesc(time,[],ks')
colorbar
colormap parula

subplot(122)
imagesc(time,[],ks'>thrks)
% colormap gray