load(fullfile(data_dir,'eeg_face_ks_ground'))

Ntrlplt = 100;
% idx = randperm(Ntrl, Ntrlplt);
load(pltidx100)
%% stat calcs
thseeg = rspdat(idx,:,:);
thsstim = stim(idx);
qstm = int16(thsstim);
crsp = copnorm_slice_omp_c_double(thseeg(:,:),Nthread);
Icop = reshape(info_c1d_slice_nobc_omp(crsp, qstm+1, 2, Ntrlplt, Nthread),[Nt Nchan]);
tt = abs(reshape(fastt2(thseeg, thsstim),[Nt Nchan]));

% permutations
Icopperm = zeros(Nt,Nchan,Nperm);
Itperm = zeros(Nt,Nchan,Nperm);
for pi=1:Nperm
    pstim = thsstim(randperm(Ntrlplt));
    qstm = int16(pstim);
    Icopperm(:,:,pi) = reshape(info_c1d_slice_nobc_omp(crsp, qstm+1, 2, Ntrlplt, Nthread),[Nt Nchan]);
    ttperm(:,:,pi) = abs(reshape(fastt2(thseeg, pstim),[Nt Nchan]));
end


%%
figure

subplot(3,4,1)
imagesc(time,[],ks')
colorbar
cl = caxis;
vline(140,'k')
vline(200,'k')
title('ks-test, 1000 trials')

subplot(3,4,2)
gtdat = 1-repmat(gt',[1 1 3]);
imagesc(time,[],gtdat)
cbax = colorbar;
set(cbax,'Visible','off')
title('permutation significance (ground truth)')

t = 140;
[~,tidx] = min(abs(time-t));
subplot(3,4,3)
topoplot(ks(tidx,:),dat.chanlocs, 'maplimits', cl);
title('140ms')

t = 200;
[~,tidx] = min(abs(time-t));
subplot(3,4,4)
topoplot(ks(tidx,:),dat.chanlocs, 'maplimits', cl);
title('200ms')


subplot(3,4,5)
imagesc(time,[],Icop')
colorbar
cl = caxis;
vline(140,'k')
vline(200,'k')
title('copula MI, 100 trials')

subplot(3,4,6)
thresh = prctile(squeeze(max(max(Icopperm,[],1),[],2)),99);
imagesc(time,[],gtdat);
hold on
imagesc(time,[],1-repmat((Icop>thresh)',[1 1 3]),'AlphaData',0.8)
cbax = colorbar;
set(cbax,'Visible','off')
title('permutation significance')

t = 140;
[~,tidx] = min(abs(time-t));
subplot(3,4,7)
topoplot(Icop(tidx,:),dat.chanlocs, 'maplimits', cl);
title('140ms')

t = 200;
[~,tidx] = min(abs(time-t));
subplot(3,4,8)
topoplot(Icop(tidx,:),dat.chanlocs, 'maplimits', cl);
title('200ms')



subplot(3,4,9)
imagesc(time,[],tt')
colorbar
cl = caxis;
vline(140,'k')
vline(200,'k')
title('t-test, 100 trials')

subplot(3,4,10)
thresh = prctile(squeeze(max(max(ttperm,[],1),[],2)),99);
imagesc(time,[],gtdat);
hold on
imagesc(time,[],1-repmat((tt>thresh)',[1 1 3]),'AlphaData',0.8)
cbax = colorbar;
set(cbax,'Visible','off')
title('permutation significance')

t = 140;
[~,tidx] = min(abs(time-t));
subplot(3,4,11)
topoplot(tt(tidx,:),dat.chanlocs, 'maplimits', cl);
title('140ms')

t = 200;
[~,tidx] = min(abs(time-t));
subplot(3,4,12)
topoplot(tt(tidx,:),dat.chanlocs, 'maplimits', cl);
title('200ms')

colormap parula