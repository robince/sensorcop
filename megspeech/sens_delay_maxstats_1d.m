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
meg = sqrt(plndat(1:Nch,:).^2 + plndat((Nch+1):end,:).^2);
meg = meg';
%% groundtruth - Spearman
tic
delays = 1:17;
Ndel = length(delays);
sp = zeros(Nch,Ndel);

tic
for di=1:Ndel
    d = delays(di);
    sp(:,di) = corr(fltspc(1:(end-d)), meg((1+d):end,:), 'type', 'spearman');
end
toc

%%
blocklen = 10 * Fs; % 10 s blocks
Nblock = ceil(Nt / blocklen);

tic
Nperm = 1000;
spperm = zeros(Nch,Ndel,Nperm);
parfor pi=1:Nperm
    thsperm = randperm(Nblock);
    for di=1:Ndel
        % align
        d = delays(di);
%         thsrsp = cmeg((1+d):end,:);
        thsspc = fltspc(1:(end-d));
        permspc = cell(1,Nblock);
        permmeg = cell(1,Nblock);
        for bi=1:Nblock
            megidx = block_index(bi,blocklen,Nt-d);
            spcidx = thsspc(block_index(thsperm(bi),blocklen,Nt-d));
            if size(megidx,2) < size(spcidx,1)1
                spcidx = spcidx(1:size(megidx,2));
            elseif size(spcidx,1) < size(megidx,2)
                megidx = megidx(:,1:size(spcidx,1));
            end
            permspc{bi} = spcidx;
            permmeg{bi} = megidx;
        end
        thsspc = cell2mat(permspc');
        megidx = cell2mat(permmeg);
        thsmeg = meg(megidx+d,:);
        spperm(:,di,pi) = corr(thsspc, thsmeg, 'type', 'spearman');
    end
end
toc

%%
maxsp = squeeze(max(max(abs(spperm),[],1),[],2));
thrsp = prctile(maxsp, 95);
gt = abs(sp)>thrsp;

save(fullfile(data_dir,'meg_speech_sp_ground'), 'sp', 'spperm', 'thrsp', 'gt','delays')

%% groundtruth - cop
tic
delays = 0:17;
Ndel = length(delays);
Icop = zeros(Nch,Ndel);

cspc = copnorm(fltspc);
cmeg = copnorm(meg);
Nthread = 16;

tic
for di=1:Ndel
    d = delays(di);
    dmeg = cmeg((1+d):end,:);
    dspc = cspc(1:(end-d));
    Icop(:,di) = info_cc_slice_nobc_omp(reshape(dmeg,[size(dmeg,1) 1 size(dmeg,2)]),1, dspc, size(dmeg,1),Nthread);
end
toc

%%
figure
subplot(121)
imagesc(sp);colorbar
subplot(122)
imagesc(Icop);colorbar

%%
blocklen = 10 * Fs; % 10 s blocks
Nblock = ceil(Nt / blocklen);

tic
Nperm = 1000;
Icopperm = zeros(Nch,Ndel,Nperm);
Nthread = 2;
parfor pi=1:Nperm
    thsperm = randperm(Nblock);
    for di=1:Ndel
        % align
        d = delays(di);
%         thsrsp = cmeg((1+d):end,:);
        thsspc = cspc(1:(end-d));
        permspc = cell(1,Nblock);
        permmeg = cell(1,Nblock);
        for bi=1:Nblock
            megidx = block_index(bi,blocklen,Nt-d);
            spcidx = thsspc(block_index(thsperm(bi),blocklen,Nt-d));
            if size(megidx,2) < size(spcidx,1)
                spcidx = spcidx(1:size(megidx,2));
            elseif size(spcidx,1) < size(megidx,2)
                megidx = megidx(:,1:size(spcidx,1));
            end
            permspc{bi} = spcidx;
            permmeg{bi} = megidx;
        end
        thsspc = cell2mat(permspc');
        megidx = cell2mat(permmeg);
        thsmeg = cmeg(megidx+d,:);
%         thsspc = cspc(spcidx);
        thsmeg = reshape(thsmeg, [size(thsmeg,1) 1 size(thsmeg,2)]);
        Icopperm(:,di,pi) =  info_cc_slice_nobc_omp(thsmeg,1, thsspc, size(thsmeg,1),Nthread);
    end
end
toc
%%
maxIcop = squeeze(max(max(Icopperm,[],1),[],2));
thrIcop = prctile(maxIcop, 95);
gt = Icop>thrIcop;

save(fullfile(data_dir,'meg_speech_Icop_ground'), 'sp', 'spperm', 'thrsp', 'gt', 'delays')

%%
gt = gtcop | gtsp;
save(fullfile(data_dir,'meg_speech_Icopsp_ground'), 'gt', 'delays');

%%
load(fullfile(data_dir,'meg_speech_Icopsp_ground'),'gt')
sampsize = [5 10 15 20 25 30]; % in blocks
Nsamps = length(sampsize);

sampres = cell(1,Nsamps);
Nperm = 100;
Nrep = 30;

blocklen = 10 * Fs; % 10 s blocks
Nblock = ceil(Nt / blocklen);


Nthread = 8;
% stats = {'Ib2' 'Ib4' 'Ib8' 'Icop' 'sp' 'pe' 'ke'};
stats = {'Ib2' 'Ib4' 'Ib8' 'Icop' 'sp' 'pe'};
Nstat = length(stats);
% cm = [];
for ii=1:Nstat
    cm.(stats{ii}) = zeros(2,2,Nrep,Nsamps);
    sig.(stats{ii}) = zeros(Nch, Ndel, Nrep, Nsamps);
end

Nthread = 4;

res = cell(Nrep,Nsamps);
for si=1:Nsamps
    thsNblock = sampsize(si);
    disp(['Nsamp: ' num2str(si)])
    tic
    repres = cell(1,Nrep);
    parfor ri=1:Nrep
        
%         if ~mod(ri,10)
%             disp(['Rep ... ' num2str(ri)])
%         end
        fres = [];
        fres.Icop = zeros(Nch,Ndel);
        fres.Ib2 = zeros(Nch,Ndel);
        fres.Ib4 = zeros(Nch,Ndel);
        fres.Ib8 = zeros(Nch,Ndel);
        fres.sp = zeros(Nch,Ndel);
        fres.pe = zeros(Nch,Ndel);
        fres.ke = zeros(Nch,Ndel);
        
        fres.Icopperm = zeros(Nch,Ndel,Nperm);
        fres.Ib2perm = zeros(Nch,Ndel,Nperm);
        fres.Ib4perm = zeros(Nch,Ndel,Nperm);
        fres.Ib8perm = zeros(Nch,Ndel,Nperm);
        fres.spperm = zeros(Nch,Ndel,Nperm);
        fres.peperm = zeros(Nch,Ndel,Nperm);
        fres.keperm = zeros(Nch,Ndel,Nperm);
        
        blkidx = randperm(Nblock, thsNblock);
        
        % quantise blocks used in this rep
        bmeg = cell(1,thsNblock);
        bspc = cell(1,thsNblock);
        blen = zeros(thsNblock,1);
        for bi=1:thsNblock
            idx = block_index(blkidx(bi),blocklen,Nt);
            bmeg{bi} = meg(idx,:);
            bspc{bi} = fltspc(idx);
            blen(bi) = length(idx);
        end
        
        thsmeg = cell2mat(bmeg');
        thsspc = cell2mat(bspc');
        
        cmeg = copnorm_slice_omp_c_double(thsmeg, Nthread);
        qmeg2 = int16(bin.eqpop_slice_omp(thsmeg,2,Nthread));
        qmeg4 = int16(bin.eqpop_slice_omp(thsmeg,4,Nthread));
        qmeg8 = int16(bin.eqpop_slice_omp(thsmeg,8,Nthread));
        
        cspc = copnorm(thsspc);
        qspc2 = int16(bin.eqpop(thsspc,2));
        qspc4 = int16(bin.eqpop(thsspc,4));
        qspc8 = int16(bin.eqpop(thsspc,8));
        
        bcmeg = cell(1,thsNblock);
        bqmeg2 = cell(1,thsNblock);
        bqmeg4 = cell(1,thsNblock);
        bqmeg8 = cell(1,thsNblock);
        bcspc = cell(1,thsNblock);
        bqspc2 = cell(1,thsNblock);
        bqspc4 = cell(1,thsNblock);
        bqspc8 = cell(1,thsNblock);
        ii = 1;
        for bi=1:thsNblock
            idx = ii:(ii+blen(bi)-1);
            ii = ii+blen(bi);
            bcmeg{bi} = cmeg(idx,:);
            bqmeg2{bi} = qmeg2(idx,:);
            bqmeg4{bi} = qmeg4(idx,:);
            bqmeg8{bi} = qmeg8(idx,:);
            
            bcspc{bi} = cspc(idx);
            bqspc2{bi} = qspc2(idx);
            bqspc4{bi} = qspc4(idx);
            bqspc8{bi} = qspc8(idx);
        end

        % TRUE VALUE
        for di=1:Ndel
            % align
            d = delays(di);

            [dmeg, dspc] = block_delay(bcmeg, bcspc, d);
            fres.Icop(:,di) = info_cc_slice_nobc_omp(reshape(dmeg,[size(dmeg,1) 1 size(dmeg,2)]),1, dspc, size(dmeg,1),Nthread);
            
            [dqmeg, dqspc] = block_delay(bqmeg2, bqspc2, d);
            fres.Ib2(:,di) = info.calc_info_slice_omp_integer_c_int16_t(dqmeg,2,dqspc,2,size(dqmeg,1),Nthread);
            [dqmeg, dqspc] = block_delay(bqmeg4, bqspc4, d);
            fres.Ib4(:,di) = info.calc_info_slice_omp_integer_c_int16_t(dqmeg,4,dqspc,4,size(dqmeg,1),Nthread);
            [dqmeg, dqspc] = block_delay(bqmeg8, bqspc8, d);
            fres.Ib8(:,di) = info.calc_info_slice_omp_integer_c_int16_t(dqmeg,8,dqspc,8,size(dqmeg,1),Nthread);
            
            [dmeg, dspc] = block_delay(bmeg, bspc, d);
            fres.sp(:,di) = abs(corr(dspc, dmeg, 'type', 'spearman'));
            fres.pe(:,di) = abs(corr(dspc, dmeg, 'type', 'pearson'));
%             fres.ke(:,di) = corr(dspc, dmeg, 'type', 'kendall');
        end
        
        % PERMS
        for pi=1:Nperm
            thsperm = randperm(thsNblock);
            for di=1:Ndel
                d = delays(di);

                [dmeg, dspc] = block_delay(bcmeg, bcspc(thsperm), d);
                fres.Icopperm(:,di,pi) = info_cc_slice_nobc_omp(reshape(dmeg,[size(dmeg,1) 1 size(dmeg,2)]),1, dspc, size(dmeg,1),Nthread);
                
                [dqmeg, dqspc] = block_delay(bqmeg2, bqspc2(thsperm), d);
                fres.Ib2perm(:,di,pi) = info.calc_info_slice_omp_integer_c_int16_t(dqmeg,2,dqspc,2,size(dqmeg,1),Nthread);
                [dqmeg, dqspc] = block_delay(bqmeg4, bqspc4(thsperm), d);
                fres.Ib4perm(:,di,pi) = info.calc_info_slice_omp_integer_c_int16_t(dqmeg,4,dqspc,4,size(dqmeg,1),Nthread);
                [dqmeg, dqspc] = block_delay(bqmeg8, bqspc8(thsperm), d);
                fres.Ib8perm(:,di,pi) = info.calc_info_slice_omp_integer_c_int16_t(dqmeg,8,dqspc,8,size(dqmeg,1),Nthread);
                
                [dmeg, dspc] = block_delay(bmeg, bspc(thsperm), d);
                fres.spperm(:,di,pi) = abs(corr(dspc, dmeg, 'type', 'spearman'));
                fres.peperm(:,di,pi) = abs(corr(dspc, dmeg, 'type', 'pearson'));
%                 fres.keperm(:,di,pi) = corr(dspc, dmeg, 'type', 'kendall');
            end
        end
        

        for ii=1:Nstat
            permmax = squeeze(max(max(fres.([stats{ii} 'perm']),[],1),[],2));
            thresh = prctile(permmax, 95);
            fres.([stats{ii} 'thresh']) = thresh;
            thssig = fres.(stats{ii}) > thresh;
            repres{ri}.(stats{ii}) = [];
            repres{ri}.(stats{ii}).sig = thssig;
            thscm = zeros(2,2);
            tp = sum(sum(thssig & gt));
            fp = sum(sum(thssig & ~gt));
            fn = sum(sum(~thssig & gt));
            tn = sum(sum(~thssig & ~gt));
            thscm(1,1) = tp;
            thscm(1,2) = fp;
            thscm(2,1) = fn;
            thscm(2,2) = tn;
            repres{ri}.(stats{ii}).cm = thscm;
        end
%         fullres{ri,si} = fres;
    end
    res(:,si) = repres;
    toc
end

%%
save res_megspeech_1d_run1 res stats sampsize