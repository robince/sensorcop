load(fullfile(data_dir,'sub2_dat.mat'));
Fs = 50;
%%
[Npln Nt] = size(plndat);

% pln sum
Nch = Npln / 2;
if round(Nch) ~= Nch
    error('Odd number of planar channels')
end
% planar sum square
megsum = sqrt(plndat(1:Nch,:).^2 + plndat((Nch+1):end,:).^2);
megsum = megsum';

megfull = zeros(Nch,2,Nt);
megfull(:,1,:) = plndat(1:Nch,:);
megfull(:,2,:) = plndat((Nch+1):end,:);
megfull = permute(megfull,[3 2 1]);

megphs = bsxfun(@rdivide, megfull, reshape(megsum,[Nt 1 Nch]));

%%
tic
delays = 0:17;
Ndel = length(delays);
Ifull = zeros(Nch,Ndel);
Isum = zeros(Nch,Ndel);
Iphs = zeros(Nch,Ndel);

cspc = copnorm(fltspc);
cmegfull = copnorm(megfull);
cmegsum = copnorm(megsum);
cmegphs = copnorm(megphs);

Nthread = 16;

tic
for di=1:Ndel
    d = delays(di);
    dspc = cspc(1:(end-d));
    
    dmegfull = cmegfull((1+d):end,:,:);
    Ifull(:,di) = info_cc_slice_nobc_omp(dmegfull, 2, dspc, size(dmegfull,1),Nthread);
    dmegsum = cmegsum((1+d):end,:);
    Isum(:,di) = info_cc_slice_nobc_omp(reshape(dmegsum,[size(dmegsum,1) 1 size(dmegsum,2)]),1, dspc, size(dmegsum,1),Nthread);
    dmegphs = cmegphs((1+d):end,:,:);
    Iphs(:,di) = info_cc_slice_nobc_omp(dmegphs, 2, dspc, size(dmegphs,1),Nthread);
end
toc


%%

blocklen = 10 * Fs; % 10 s blocks
Nblock = ceil(Nt / blocklen);

tic
Nperm = 200;
Ifullperm = zeros(Nch,Ndel,Nperm);
Iphsperm = zeros(Nch,Ndel,Nperm);
Isumperm = zeros(Nch,Ndel,Nperm);

% build blocks
bmegfull = cell(1,Nblock);
bmegspc = cell(1,Nblock);
bmegphs = cell(1,Nblock);
bspc = cell(1,Nblock);
blen = zeros(Nblock,1);
for bi=1:Nblock
    idx = block_index(bi,blocklen,Nt);
    bmegfull{bi} = cmegfull(idx,:,:);
    bmegphs{bi} = cmegphs(idx,:,:);
    bmegsum{bi} = cmegsum(idx,:);
    bspc{bi} = cspc(idx);
    blen(bi) = length(idx);
end

Nthread = 4;
parfor pi=1:Nperm
    thsperm = randperm(Nblock);
    for di=1:Ndel
        d = delays(di);
        
        [dmegsum, dspc] = block_delay(bmegsum, bspc(thsperm), d);
        thsmeg = reshape(dmegsum, [size(dmegsum,1) 1 size(dmegsum,2)]);
        Isumperm(:,di,pi) =  info_cc_slice_nobc_omp(thsmeg,1, dspc, size(thsmeg,1),Nthread);
        
        [dmegfull, dspc] = block_delay2d(bmegfull, bspc(thsperm), d);
        Ifullperm(:,di,pi) =  info_cc_slice_nobc_omp(dmegfull, 2, dspc, size(dmegfull,1),Nthread);

        [dmegphs, dspc] = block_delay2d(bmegphs, bspc(thsperm), d);
        Iphsperm(:,di,pi) =  info_cc_slice_nobc_omp(dmegphs, 2, dspc, size(dmegphs,1),Nthread);

    end
end
toc
%%
thrprc = 99;
maxIfull = squeeze(max(max(Ifullperm,[],1),[],2));
thrIfull = prctile(maxIfull,thrprc);

maxIsum = squeeze(max(max(Isumperm,[],1),[],2));
thrIsum = prctile(maxIsum,thrprc);

maxIphs = squeeze(max(max(Iphsperm,[],1),[],2));
thrIphs = prctile(maxIphs,thrprc);

%%
figure
subplot(3,3,1)
di = 8;
imagesc(1000*(delays/Fs),[],Ifull);
vline(1000*delays(di)/Fs,'k');
colorbar
cl = caxis;
title('2d planar gradient')
subplot(3,3,2)
imagesc(1000*(delays/Fs),[],1-repmat(Ifull>thrIfull,[1 1 3]))
title('permutation significance')
subplot(3,3,3)
topoplot(Ifull(:,di)',chanlocs,'maplimits',cl);
title(sprintf('%d ms', 1000*delays(di)/Fs));

subplot(3,3,4)
di = 6;
imagesc(1000*(delays/Fs),[],Isum);
vline(1000*delays(di)/Fs,'k');
colorbar
cl = caxis;
title('planar gradient amplitude')
subplot(3,3,5)
imagesc(1000*(delays/Fs),[],1-repmat(Isum>thrIsum,[1 1 3]))
title('permutation significance')
subplot(3,3,6)
topoplot(Isum(:,di)',chanlocs,'maplimits',cl);
title(sprintf('%d ms', 1000*delays(di)/Fs));

subplot(3,3,7)
di = 8;
imagesc(1000*(delays/Fs),[],Iphs)
vline(1000*delays(di)/Fs,'k');
colorbar
cl = caxis;
title('planar gradient direction')
subplot(3,3,8)
imagesc(1000*(delays/Fs),[],1-repmat(Iphs>thrIphs,[1 1 3]))
title('permutation significance')
subplot(3,3,9)
topoplot(Iphs(:,di)',chanlocs,'maplimits',cl);
title(sprintf('%d ms', 1000*delays(di)/Fs));

colormap parula