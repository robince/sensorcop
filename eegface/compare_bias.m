% compare bias between different info methods

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
% chi = 41;
% ti = 67;

chi = 14;
ti = 45;


log_samp_size = 3:10;
Nsamp = length(log_samp_size);

Nrep = 600;
res.Ib2 = zeros(Nrep, Nsamp);
res.Ib4 = zeros(Nrep, Nsamp);
res.Ib8 = zeros(Nrep, Nsamp);
res.Ib2mm = zeros(Nrep, Nsamp);
res.Ib4mm = zeros(Nrep, Nsamp);
res.Ib8mm = zeros(Nrep, Nsamp);
res.Icop = zeros(Nrep, Nsamp);
Nthread = 4;
for sampi=1:Nsamp
    thsNtrl = 2.^(log_samp_size(sampi));
    for repi=1:Nrep
        idx = randperm(Ntrl, thsNtrl);
        thseeg = rspdat(idx,ti,chi);
        thsstim = stim(idx);
        qstm = int16(thsstim);
        
        % true values
%         qrsp2 = int16(bin.eqpop_slice_omp(thseeg(:,:),2,Nthread));
        qrsp2 = int16(bin.eqpop(thseeg,2));
        res.Ib2(repi,sampi) = info.calc_info_slice_omp_integer_c_int16_t(qrsp2, 2, qstm, 2, thsNtrl, Nthread);
        res.Ib2mm(repi,sampi) = info.calc_info_slice_omp_integer_c_int16_t(qrsp2, 2, qstm, 2, thsNtrl, Nthread) - 1./(2*thsNtrl*log(2));

%         qrsp4 = int16(bin.eqpop_slice_omp(thseeg(:,:),4,Nthread));
        qrsp4 = int16(bin.eqpop_slice(thseeg,4));
        res.Ib4(repi,sampi) = info.calc_info_slice_omp_integer_c_int16_t(qrsp4, 4, qstm, 2, thsNtrl, Nthread);
        res.Ib4mm(repi,sampi) = info.calc_info_slice_omp_integer_c_int16_t(qrsp4, 4, qstm, 2, thsNtrl, Nthread) - 3./(2*thsNtrl*log(2));

%         qrsp8 = int16(bin.eqpop_slice_omp(thseeg(:,:),8,Nthread));
        qrsp8 = int16(bin.eqpop_slice(thseeg,8));
        res.Ib8(repi,sampi) = info.calc_info_slice_omp_integer_c_int16_t(qrsp8, 8, qstm, 2, thsNtrl, Nthread);
        res.Ib8mm(repi,sampi) = info.calc_info_slice_omp_integer_c_int16_t(qrsp8, 8, qstm, 2, thsNtrl, Nthread) - 7./(2*thsNtrl*log(2));
        
        crsp = copnorm(thseeg);
        try
            res.Icop(repi,sampi) = info_gd(crsp, qstm, 2, true, false, false);
        catch
            res.Icop(repi,sampi) = NaN;
        end
    end
end

%%
est = {'Ib2' 'Ib4' 'Ib8' 'Ib2mm' 'Ib4mm' 'Ib8mm' 'Icop'};
figure
hold all
for ei=1:length(est)
    plot(log_samp_size, nanmean(res.(est{ei}),1));
end
legend(est)
%%
est = {'Ib2' 'Ib4' 'Ib8' 'Icop'};
figure
hold all
for ei=1:length(est)
    plot(log_samp_size, nanstd(res.(est{ei}),1));
end
legend(est)
%%
est = {'Ib4'  'Icop'};
figure
hold all
for ei=1:length(est)
    errorbar(2.^log_samp_size, nanmean(res.(est{ei}),1),nanstd(res.(est{ei}),[],1));
end
legend(est)