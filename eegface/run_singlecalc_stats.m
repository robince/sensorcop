% run trial-dependendent permutation tests of various statistical tests

subid = 'HKI1';
fname = sprintf('%s_ridat.mat',subid);
dat = load(fullfile(data_dir, fname));

time = dat.time;
stim = dat.stim;
% rspdat = permute(dat.eegdat, [3 2 1]);
rspdat = permute(dat.csddat, [3 2 1]);
[Ntrl Nt Nchan] = size(rspdat);

%%
% max effect
ci = 41; % B9
ti = 319;
% ti = 294;

fname = sprintf('%s_csd_%d_%d_trialselect.mat',subid,ci,ti);
load(fullfile(data_dir,fname));
%%
rsp = rspdat(:,ti,ci);
rsp0 = rsp(stim==0);
rsp1 = rsp(stim==1);
trlnums = [8 16 32 64 128 256 512 1024];
Ntrlnum = length(trlnums);
Nperm = 500;

% select trials
Ntrlrep = 100;
trldat = cell(Ntrlrep,Ntrlnum);
trlstm = cell(1,Ntrlnum);
trlperm = cell(1,Ntrlnum);
for i=1:Ntrlnum
    Nts = trlnums(i) / 2;
    for ri=1:Ntrlrep
%         trldat{ri,i} = [ rsp0(randperm(length(rsp0),Nts)); rsp1(randperm(length(rsp1),Nts)) ];        
        trldat{ri,i} = [ rsp0(randi(length(rsp0),[Nts 1])); rsp1(randi(length(rsp1),[Nts 1])) ];
    end
    trlstm{i} = [ zeros(Nts,1); ones(Nts,1) ];
    
    permdat = zeros(trlnums(i),Nperm);
    for pi=1:Nperm
        % without replacement
        permdat(:,pi) = randperm(trlnums(i));
        % with replacement
%          permdat(:,pi) = randi(trlnums(i),[trlnums(i) 1]);
    end
    trlperm{i} = permdat;
end

fname = sprintf('%s_csd_%d_%d_trialselect.mat',subid,ci,ti);
save(fullfile(data_dir,fname), 'trldat','trlstm','trlperm')

%% 1D RESPONSE
% ci = 41; % B9
% % ti = 319;
% ti = 294;
fname = sprintf('%s_csd_%d_%d_trialselect.mat',subid,ci,ti);
load(fullfile(data_dir,fname), 'trldat','trlstm','trlperm')

Ntrlnum = length(trlstm);
Nperm = size(trlperm{1},2);
Ntrlrep = size(trldat,1);

tic
trlres = cell(1,Ntrlnum);
parfor i=1:Ntrlnum
    % quantise
    res = [];
    res.Ic = zeros(1,Ntrlrep);
    res.Icperm = zeros(Nperm,Ntrlrep);
    res.Ig = zeros(1,Ntrlrep);
    res.Igperm = zeros(Nperm,Ntrlrep);
    res.Ib2 = zeros(1,Ntrlrep);
    res.Ib2perm = zeros(Nperm,Ntrlrep);
    res.Ib4 = zeros(1,Ntrlrep);
    res.Ib4perm = zeros(Nperm,Ntrlrep);
    res.Ib8 = zeros(1,Ntrlrep);
    res.Ib8perm = zeros(Nperm,Ntrlrep);
    res.t = zeros(1,Ntrlrep);
    res.tperm = zeros(Nperm,Ntrlrep);
    res.t2 = zeros(1,Ntrlrep);
    res.t2perm = zeros(Nperm,Ntrlrep);
    res.ks = zeros(1,Ntrlrep);
    res.ksperm = zeros(Nperm,Ntrlrep);
%     res.logregaz = zeros(1,Ntrlrep);
%     res.logregazperm = zeros(1,Ntrlrep);
    
    Ntrl = size(trldat{1,i},1);
    for ri=1:Ntrlrep
        rdat = trldat{ri,i};
        cdat = copnorm(rdat);
        qdat2 = int16(bin.eqpop(rdat,2));
        qdat4 = int16(bin.eqpop(rdat,4));
        qdat8 = int16(bin.eqpop(rdat,8));
        
        % 2 sample format
        rdat0 = rdat(1:(Ntrl/2));
        rdat1 = rdat(((Ntrl/2)+1):end);
        
        % true val
        res.Ic(ri) = info_gd(cdat, trlstm{i}, 2, true, true, false);
        res.Ig(ri) = info_gd(rdat, trlstm{i}, 2, true, false, false);
        qstm = int16(trlstm{i});
        res.Ib2(ri) = info.calc_info_integer_c_int16_t(qdat2,2,qstm,2,Ntrl);
        res.Ib4(ri) = info.calc_info_integer_c_int16_t(qdat4,4,qstm,2,Ntrl);
        res.Ib8(ri) = info.calc_info_integer_c_int16_t(qdat8,8,qstm,2,Ntrl);
        
        [h, p, ci, stat] = ttest2(rdat0,rdat1);
        res.t(ri) = stat.tstat;
        [h, p, ksstat] = kstest2(rdat0,rdat1);
        res.ks(ri) = ksstat;
        res.t2(ri) = res.t(ri).^2;
        
%         model = logregFit(rdat, trlstm{i});
%         [yhat, p] = logregPredict(model, rdat);
%         res.logregaz(ri) = rocarea_mex(p, int16(trlstm{i}));

        % permutations
        for pi=1:Nperm
            stm = trlstm{i}(trlperm{i}(:,pi));
            rdat0 = rdat(stm==0);
            rdat1 = rdat(stm==1);
            qstm = int16(stm);
            res.Icperm(pi,ri) = info_gd(cdat, stm, 2, true, true, false);
            res.Igperm(pi,ri) = info_gd(rdat, stm, 2, true, false, false);
            res.Ib2perm(pi,ri) = info.calc_info_integer_c_int16_t(qdat2,2,qstm,2,Ntrl);
            res.Ib4perm(pi,ri) = info.calc_info_integer_c_int16_t(qdat4,4,qstm,2,Ntrl);
            res.Ib8perm(pi,ri) = info.calc_info_integer_c_int16_t(qdat8,8,qstm,2,Ntrl);
            [h, p, ci, stat] = ttest2(rdat0,rdat1);
            res.tperm(pi,ri) = stat.tstat;
            [h, p, ksstat] = kstest2(rdat0,rdat1);
            res.ksperm(pi,ri) = ksstat;
            res.t2perm(pi,ri) = res.tperm(pi,ri).^2;
        end
    end
    trlres{i} = res;
end
toc
res1d = trlres;

%% 2D RESPONSE
% ci = 41; % B9
% % ti = 319;
% ti = 294;
fname = sprintf('%s_csd_%d_%d_trialselect.mat',subid,ci,ti);
load(fullfile(data_dir,fname), 'trldat','trlstm','trlperm')


%%
rrsp = rspdat(:,ti,ci);
drsp = rspdat(:,ti+1,ci) - rspdat(:,ti-1,ci);
rsp = [rrsp drsp];

rsp0 = rsp(stim==0,:);
rsp1 = rsp(stim==1,:);
trlnums = [8 16 32 64 128 256 512 1024];
Ntrlnum = length(trlnums);


% select trials
Ntrlrep = 100;
trldat2d = cell(Ntrlrep,Ntrlnum);

for i=1:Ntrlnum
    Nts = trlnums(i) / 2;
    for ri=1:Ntrlrep
%         trldat{ri,i} = [ rsp0(randperm(length(rsp0),Nts)); rsp1(randperm(length(rsp1),Nts)) ];        
        trldat2d{ri,i} = [ rsp0(randi(length(rsp0),[Nts 1]),:); rsp1(randi(length(rsp1),[Nts 1]),:) ];
    end
    trlstm{i} = [ zeros(Nts,1); ones(Nts,1) ];
    
    permdat = zeros(trlnums(i),Nperm);
    for pi=1:Nperm
        % without replacement
        permdat(:,pi) = randperm(trlnums(i));
        % with replacement
%          permdat(:,pi) = randi(trlnums(i),[trlnums(i) 1]);
    end
    trlperm{i} = permdat;
end

fname = sprintf('%s_csd_%d_%d_trialselect.mat',subid,ci,ti);
save('-append',fullfile(data_dir,fname), 'trldat2d')

%%
Ntrlnum = length(trlstm);
Nperm = size(trlperm{1},2);
Ntrlrep = size(trldat,1);

tic
trlres = cell(1,Ntrlnum);
parfor i=1:Ntrlnum
    % quantise
    res = [];
    res.Ic = zeros(1,Ntrlrep);
    res.Icperm = zeros(Nperm,Ntrlrep);
    res.Ig = zeros(1,Ntrlrep);
    res.Igperm = zeros(Nperm,Ntrlrep);
    res.Ib2 = zeros(1,Ntrlrep);
    res.Ib2perm = zeros(Nperm,Ntrlrep);
    res.Ib4 = zeros(1,Ntrlrep);
    res.Ib4perm = zeros(Nperm,Ntrlrep);
    res.Ib8 = zeros(1,Ntrlrep);
    res.Ib8perm = zeros(Nperm,Ntrlrep);
    res.ht2 = zeros(1,Ntrlrep);
    res.ht2perm = zeros(Nperm,Ntrlrep);
%     res.logregaz = zeros(1,Ntrlrep);
%     res.logregazperm = zeros(1,Ntrlrep);
    
    Ntrl = size(trldat{1,i},1);
    for ri=1:Ntrlrep
        rdat = trldat2d{ri,i};
        cdat = copnorm(rdat);
        qdat2 = int16(numbase2dec(double(bin.eqpop_slice(rdat,2))',2))';
        qdat4 = int16(numbase2dec(double(bin.eqpop_slice(rdat,4))',4))';
        qdat8 = int16(numbase2dec(double(bin.eqpop_slice(rdat,8))',8))';
        
        % 2 sample format
        rdat0 = rdat(1:(Ntrl/2),:);
        rdat1 = rdat(((Ntrl/2)+1):end,:);
        
        % true val
        res.Ic(ri) = info_gd(cdat, trlstm{i}, 2, true, true, false);
        res.Ig(ri) = info_gd(rdat, trlstm{i}, 2, true, false, false);
        qstm = int16(trlstm{i});
        res.Ib2(ri) = info.calc_info_integer_c_int16_t(qdat2,4,qstm,2,Ntrl);
        res.Ib4(ri) = info.calc_info_integer_c_int16_t(qdat4,16,qstm,2,Ntrl);
        res.Ib8(ri) = info.calc_info_integer_c_int16_t(qdat8,64,qstm,2,Ntrl);
        
%         [h, p, ci, stat] = ttest2(rdat0,rdat1);
%         res.t(ri) = stat.tstat;
%         [h, p, ksstat] = kstest2(rdat0,rdat1);
%         res.ks(ri) = ksstat;
        res.ht2(ri) = ht2test(rdat0,rdat1);
%         model = logregFit(rdat, trlstm{i});
%         [yhat, p] = logregPredict(model, rdat);
%         res.logregaz(ri) = rocarea_mex(p, int16(trlstm{i}));

        % permutations
        for pi=1:Nperm
            stm = trlstm{i}(trlperm{i}(:,pi));
            rdat0 = rdat(stm==0,:);
            rdat1 = rdat(stm==1,:);
            qstm = int16(stm);
            res.Icperm(pi,ri) = info_gd(cdat, stm, 2, true, true, false);
            res.Igperm(pi,ri) = info_gd(rdat, stm, 2, true, false, false);
            res.Ib2perm(pi,ri) = info.calc_info_integer_c_int16_t(qdat2,2,qstm,2,Ntrl);
            res.Ib4perm(pi,ri) = info.calc_info_integer_c_int16_t(qdat4,16,qstm,2,Ntrl);
            res.Ib8perm(pi,ri) = info.calc_info_integer_c_int16_t(qdat8,64,qstm,2,Ntrl);
%             [h, p, ci, stat] = ttest2(rdat0,rdat1);
%             res.tperm(pi,ri) = stat.tstat;
%             [h, p, ksstat] = kstest2(rdat0,rdat1);
%             res.ksperm(pi,ri) = ksstat;
            res.ht2perm(pi,ri) = ht2test(rdat0,rdat1);
        end
    end
    trlres{i} = res;
end
toc
res2d = trlres;


%% z-score stats
trlres = res1d;
stats = {'Ic' 'Ig' 'Ib2' 'Ib4' 'Ib8' 't' 'ks'};
cubetransform = [1 1 1 1 1 0 0];
% cubetransform = [0 0 0 0 0 0 0];
% % 
% trlres = res2d;
% stats = {'Ic' 'Ig' 'Ib2' 'Ib4' 'Ib8' 'ht2'};
% cubetransform = [1 1 1 1 1 1];

Z = 0; % normal z-score
% Z = 1; % median + NIQR

for i=1:Ntrlnum
    res = trlres{i};
    for sti=1:length(stats)
        stat = stats{sti};
        
        res.([stat 'Z']) = zeros(1,Ntrlrep);
        for ri=1:Ntrlrep
            if cubetransform(sti)
                perm = res.([stat 'perm'])(:,ri);
                offset = min([min(perm) 0]);
                perm = (perm-offset).^(1/3);
                val = (res.(stat)(ri)-offset).^(1/3);
            else
                % no transform
                perm = res.([stat 'perm'])(:,ri);
                val = res.(stat)(ri);
            end
                
            if Z==0
                mu = mean(perm);
                sigma = std(perm);
            elseif Z==1
                mu = median(perm);
                sigma = 0.7413*iqr(perm);
            else
                clear mu sigma
            end
            res.([stat 'Z'])(ri) = (val - mu) / sigma;
        end
    end
%     if isfield(res,'ht2perm')
%         res.htZ = zeros(1,Ntrlrep);
%         for ri=1:Ntrlrep
%             if Z==0
%                 mu = mean(sqrt(res.ht2perm(:,ri)));
%                 sigma = std(sqrt(res.ht2perm(:,ri)));
%             elseif Z==1
%                 mu = median(sqrt(res.ht2perm(:,ri)));
%                 sigma = 0.7413*iqr(sqrt(res.ht2perm(:,ri)));
%             else
%                 clear mu sigma
%             end
%             res.htZ(ri) = (sqrt(res.ht2(ri)) - mu) / sigma;
%         end
%     end
%     if isfield(res,'tperm')
%         res.t2Z = zeros(1,Ntrlrep);
%         for ri=1:Ntrlrep
%             if Z==0
%                 mu = mean(res.tperm(:,ri).^2);
%                 sigma = std(res.tperm(:,ri).^2);
%             elseif Z==1
%                 mu = median(res.tperm(:,ri).^2);
%                 sigma = 0.7413*iqr(res.tperm(:,ri).^2);
%             else
%                 clear mu sigma
%             end
%             res.t2Z(ri) = (res.t(ri).^2 - mu) / sigma;
%         end
%     end
    trlres{i} = res;
end


%% dat for graph
% stats = {'Ic' 'Ig' 'Ib2' 'Ib4' 'Ib8' 't' 'ks'};
stats = {'Ic' 'Ib2' 'Ib4' 'Ib8' 't' 'ks'};

% stats = {'Ic'  'Ib4' 'Ib2'  't' 'ks'};
% stats = {'Ic'  'Ib4' 'Ib2'  'ht' 'ht2' };
pltavg = [];
pltstd = [];

for sti=1:length(stats)
    stat = stats{sti};
    
    pltavg.(stat) = zeros(1,Ntrlnum);
    pltstd.(stat) = zeros(1,Ntrlnum);
    
    for i=1:Ntrlnum   
        pltavg.(stat)(i) = mean(abs(trlres{i}.([stat 'Z'])));
        pltstd.(stat)(i) = std(abs(trlres{i}.([stat 'Z'])));
    end
end

figure
for sti=1:length(stats)
    stat = stats{sti};
    loglog(trlnums, pltavg.(stat))
%     semilogx(trlnums, pltavg.(stat))
    hold all
end
legend(stats)
set(gca,'TickDir','out')

%% hypothesis reject?
trlres = res1d;
stats = {'Ic' 'Ig' 'Ib2' 'Ib4' 'Ib8' 't' 'ks'};

% % 
% trlres = res2d;
% stats = {'Ic' 'Ig' 'Ib2' 'Ib4' 'Ib8' 'ht2'};

for i=1:Ntrlnum
    res = trlres{i};
    for sti=1:length(stats)
        stat = stats{sti};
        
        res.([stat 'H']) = zeros(1,Ntrlrep);
        if strcmpi(stat,'t')
            % two tailed            
            for ri=1:Ntrlrep
                thr_lo = prctile( res.([stat 'perm'])(:,ri), 0.5);
                thr_high = prctile( res.([stat 'perm'])(:,ri), 99.5);
                res.([stat 'H'])(ri) = (res.(stat)(ri) > thr_high) || (res.(stat)(ri) < thr_lo);
            end
        else
            % one tailed
            for ri=1:Ntrlrep
                thr = prctile( res.([stat 'perm'])(:,ri), 99);
                res.([stat 'H'])(ri) = res.(stat)(ri) > thr;
            end
        end
    end

    trlres{i} = res;
end

%% plot H
% stats = {'Ic' 'Ig' 'Ib2' 'Ib4' 'Ib8' 't' 'ks'};
stats = {'Ic' 'Ib2' 'Ib4' 'Ib8' 't' 'ks'};

% stats = {'Ic'  'Ib4' 'Ib2'  't' 'ks'};
% stats = {'Ic'  'Ib4' 'Ib2'  'ht' 'ht2' };
pltavg = [];


for sti=1:length(stats)
    stat = stats{sti};
    
    pltavg.(stat) = zeros(1,Ntrlnum);
    
    for i=1:Ntrlnum   
        pltavg.(stat)(i) = mean(trlres{i}.([stat 'H']));
    end
end
figure
for sti=1:length(stats)
    stat = stats{sti};
%     loglog(trlnums, pltavg.(stat))
    semilogx(trlnums, pltavg.(stat))
    hold all
end
legend(stats)
set(gca,'TickDir','out')


%%
stats = {'Ic' 'Ig' 'Ib2' 'Ib4' 'Ib8'};
stats = {'Ic'  'Ib2' 'Ib4' 'Ib8'  };
pltavg = [];
pltstd = [];

for sti=1:length(stats)
    stat = stats{sti};
    
    pltavg.(stat) = zeros(1,Ntrlnum);
    pltstd.(stat) = zeros(1,Ntrlnum);
    
    for i=1:Ntrlnum   
        pltavg.(stat)(i) = mean(abs(trlres{i}.([stat])));
        pltstd.(stat)(i) = std(abs(trlres{i}.([stat])));
    end
end
figure
for sti=1:length(stats)
    stat = stats{sti};
    semilogx(trlnums, pltavg.(stat))
    hold all
end
legend(stats)
%%
stats = {'t' 'ks'};
pltavg = [];
pltstd = [];

for sti=1:length(stats)
    stat = stats{sti};
    
    pltavg.(stat) = zeros(1,Ntrlnum);
    pltstd.(stat) = zeros(1,Ntrlnum);
    
    for i=1:Ntrlnum   
        pltavg.(stat)(i) = mean(abs(trlres{i}.([stat])));
        pltstd.(stat)(i) = std(abs(trlres{i}.([stat])));
    end
end
figure
for sti=1:length(stats)
    stat = stats{sti};
    semilogx(trlnums, pltavg.(stat))
    hold all
end
legend(stats)
 
