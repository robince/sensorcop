subid = 'HKI1';
fname = sprintf('%s_ridat.mat',subid);
dat = load(fullfile(data_dir, fname));

time = dat.time;
stim = dat.stim;
% rspdat = permute(dat.eegdat, [3 2 1]);
rspdat = permute(dat.csddat, [3 2 1]);
[Ntrl Nt Nchan] = size(rspdat);

tidx = (time>0) & (time<500);
time = time(tidx);
rspdat = rspdat(:,tidx,:);
% just keep 1000 trials
rspdat = rspdat(1:1000,:,:);
[Ntrl Nt Nchan] = size(rspdat);
stim = stim(1:1000);
%% groundtruth - KS test
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
thrks = prctile(maxks, 95);
gt = ks>thrks;

save(fullfile(data_dir,'eeg_face_ks_ground'), 'ks', 'ksperm', 'thrks', 'gt')

%%
load(fullfile(data_dir,'eeg_face_ks_ground'),'gt')
sampsize = [25 50 100 200 500];
Nsamps = length(sampsize);

sampres = cell(1,Nsamps);
Nperm = 200;
Nrep = 50;

Nthread = 32;
stats = {'Ib2' 'Ib4' 'Ib8' 'Icop' 't' 'ks'};
stats = {'Ib8'}
Nstat = length(stats);
% cm = [];
for ii=1:Nstat
    cm.(stats{ii}) = zeros(2,2,Nrep,Nsamps);
end
fullres = cell(Nrep,Nsamps);

for si=1:Nsamps
%     Nfold = 1000./sampsize(si);
%     ind = crossvalind('Kfold', Ntrl, Nfold);
    thsNtrl = sampsize(si);
    disp(['Nsamp: ' num2str(si)])
    tic
    for ri=1:Nrep
%         if ~mod(ri,10)
            disp(['Rep ... ' num2str(ri)])
%         end
        fres = [];
%         idx = ind==fi;
        idx = randperm(Ntrl, thsNtrl);
        thseeg = rspdat(idx,:,:);
        thsstim = stim(idx);
        qstm = int16(thsstim);
        
        % true values
        qrsp2 = int16(bin.eqpop_slice_omp(thseeg(:,:),2,Nthread));
        fres.Ib2 = reshape(info.calc_info_slice_omp_integer_c_int16_t(qrsp2, 2, qstm, 2, thsNtrl, Nthread), [Nt Nchan]);

        qrsp4 = int16(bin.eqpop_slice_omp(thseeg(:,:),4,Nthread));
        fres.Ib4 = reshape(info.calc_info_slice_omp_integer_c_int16_t(qrsp4, 4, qstm, 2, thsNtrl, Nthread), [Nt Nchan]);

        qrsp8 = int16(bin.eqpop_slice_omp(thseeg(:,:),8,Nthread));
        fres.Ib8 = reshape(info.calc_info_slice_omp_integer_c_int16_t(qrsp8, 8, qstm, 2, thsNtrl, Nthread), [Nt Nchan]);
        
        crsp = copnorm_slice_omp_c_double(thseeg(:,:),Nthread);
        fres.Icop = reshape(info_c1d_slice_nobc_omp(crsp, qstm+1, 2, thsNtrl, Nthread),[Nt Nchan]);
        
        fres.t = abs(reshape(fastt2(thseeg, thsstim),[Nt Nchan]));
        
        fres.ks = reshape(kstest_slice_omp(thseeg(:,:), qstm, thsNtrl, Nthread), [Nt Nchan]);

        
        % permutations
        Ib2perm = zeros(Nt,Nchan,Nperm);
        Ib4perm = zeros(Nt,Nchan,Nperm);
        Ib8perm = zeros(Nt,Nchan,Nperm);
        Icopperm = zeros(Nt,Nchan,Nperm);
        Iksperm = zeros(Nt,Nchan,Nperm);
        Itperm = zeros(Nt,Nchan,Nperm);
        for pi=1:Nperm
            pstim = thsstim(randperm(thsNtrl));
            qstm = int16(pstim);
%             fres.Ib2perm(:,:,pi) = reshape(info.calc_info_slice_omp_integer_c_int16_t(qrsp2, 2, qstm, 2, thsNtrl, Nthread), [Nt Nchan]);
%             fres.Ib4perm(:,:,pi) = reshape(info.calc_info_slice_omp_integer_c_int16_t(qrsp4, 4, qstm, 2, thsNtrl, Nthread), [Nt Nchan]);
            fres.Ib8perm(:,:,pi) = reshape(info.calc_info_slice_omp_integer_c_int16_t(qrsp8, 8, qstm, 2, thsNtrl, Nthread), [Nt Nchan]);
%             fres.Icopperm(:,:,pi) = reshape(info_c1d_slice_nobc_omp(crsp, qstm+1, 2, thsNtrl, Nthread),[Nt Nchan]);
%             fres.tperm(:,:,pi) = abs(reshape(fastt2(thseeg, pstim),[Nt Nchan]));
%             fres.ksperm(:,:,pi) = reshape(kstest_slice_omp(thseeg(:,:), qstm, thsNtrl, Nthread), [Nt Nchan]);
        end

        for ii=1:Nstat
            permmax = squeeze(max(max(fres.([stats{ii} 'perm']),[],1),[],2));
            thresh = prctile(permmax, 95);
            fres.([stats{ii} 'thresh']) = thresh;
            sig = fres.(stats{ii}) > thresh;
            tp = sum(sum(sig & gt));
            fp = sum(sum(sig & ~gt));
            fn = sum(sum(~sig & gt));
            tn = sum(sum(~sig & ~gt));
            cm.(stats{ii})(1,1,ri,si) = tp;
            cm.(stats{ii})(1,2,ri,si) = fp;
            cm.(stats{ii})(2,1,ri,si) = fn;
            cm.(stats{ii})(2,2,ri,si) = tn;
        end
%         fullres{ri,si} = fres;
    end
    toc
end

    
