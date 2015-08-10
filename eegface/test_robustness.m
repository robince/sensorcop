% compare robustness of t-test vs Icop
data_dir = '~/Documents/gladata/sensorcop';
addpath('~/Documents/glacode/para_info/')
addpath('~/Documents/glacode/para_info/mex')
addpath('~/Documents/glacode/info')
%%
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
[Ntrl Nt Nchan] = size(rspdat);

%%
load(fullfile(data_dir,'eeg_face_ks_ground'),'gt')
chi = 41;
ti = 67;

thsNtrl = 100;
corrupt_prct = [0 5 10 15 20 25 30 35 40 45 50];
Ncorr = length(corrupt_prct);

Nrep = 100;
Ntrl = 100;

% calculate value with X% of outliers added (from random gaussian
% with 5*SD of data.



corrres = cell(1,Ncorr);
Nrep = 50;
Nperm = 200;

% Nthread = 32;
Nthread = 8;
% stats = {'Ib2' 'Ib4' 'Ib8' 'Icop' 't' 'ks'};
stats = {'Ib2' 'Ib4' 'Ib8' 'Icop' 't'};

Nstat = length(stats);
% cm = [];
for ii=1:Nstat
    cm.(stats{ii}) = zeros(2,2,Nrep,Ncorr);
end

for si=1:Ncorr
    disp(['Ncorr: ' num2str(si)])
    tic
    for ri=1:Nrep
        if ~mod(ri,10)
            disp(['Rep ... ' num2str(ri)])
        end
        fres = [];
        
        idx = randperm(Ntrl, thsNtrl);
        thseeg = rspdat(idx,:,:);
        thsstim = stim(idx);
        qstm = int16(thsstim);
        
        corridx = randperm(thsNtrl, round(thsNtrl*corrupt_prct(sampi)/100));
        thseeg(corridx,:,:) = 5.*bsxfun(@times, randn(length(corridx),size(thseeg,2),size(thseeg,3)), std(thseeg,[],1));

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
        
%         fres.ks = reshape(kstest_slice_omp(thseeg(:,:), qstm, thsNtrl, Nthread), [Nt Nchan]);

        
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
            fres.Ib2perm(:,:,pi) = reshape(info.calc_info_slice_omp_integer_c_int16_t(qrsp2, 2, qstm, 2, thsNtrl, Nthread), [Nt Nchan]);
            fres.Ib4perm(:,:,pi) = reshape(info.calc_info_slice_omp_integer_c_int16_t(qrsp4, 4, qstm, 2, thsNtrl, Nthread), [Nt Nchan]);
            fres.Ib8perm(:,:,pi) = reshape(info.calc_info_slice_omp_integer_c_int16_t(qrsp8, 8, qstm, 2, thsNtrl, Nthread), [Nt Nchan]);
            fres.Icopperm(:,:,pi) = reshape(info_c1d_slice_nobc_omp(crsp, qstm+1, 2, thsNtrl, Nthread),[Nt Nchan]);
            fres.tperm(:,:,pi) = abs(reshape(fastt2(thseeg, pstim),[Nt Nchan]));
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
        
