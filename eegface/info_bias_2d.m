% compare bias between different info methods
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
% add gradient
rspdat = permute(rspdat,[2 1 3]);
drsp = gradient_dim1(rspdat);
rspdat = reshape(permute(rspdat,[2 1 3]),[Ntrl 1 Nt Nchan]);
drsp = reshape(permute(drsp,[2 1 3]),[Ntrl 1 Nt Nchan]);
rspdat = cat(2, rspdat, drsp);

%%
chi = 41;
ti = 67;

% chi = 14;
% ti = 45;


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
    thsNtrl = 2.^(log_samp_size(sampi))
    for repi=1:Nrep
        idx = randperm(Ntrl, thsNtrl);
        thsstim = stim(idx);
        ss = sum(thsstim==0);
        while (ss<2 || (thsNtrl-ss)<2)
            idx = randperm(Ntrl, thsNtrl);
            thsstim = stim(idx);
            ss = sum(thsstim==0);
        end
        thseeg = rspdat(idx,:,ti,chi);
        qstm = int16(thsstim);
        
        % true values
        rqrsp2 = double(bin.eqpop(squeeze(thseeg(:,1)),2));
        dqrsp2 = double(bin.eqpop(squeeze(thseeg(:,2)),2));
        qrsp2 = numbase2dec([rqrsp2(:) dqrsp2(:)]',2);
        qrsp2 = int16(reshape(qrsp2, [thsNtrl 1]));
        res.Ib2(repi,sampi) = info.calc_info_integer_c_int16_t(qrsp2, 4, qstm, 2, thsNtrl);
        res.Ib2mm(repi,sampi) = res.Ib2(repi,sampi) - 3./(2*thsNtrl*log(2));

        rqrsp4 = double(bin.eqpop(squeeze(thseeg(:,1)),4));
        dqrsp4 = double(bin.eqpop(squeeze(thseeg(:,2)),4));
        qrsp4 = numbase2dec([rqrsp4(:) dqrsp4(:)]',4);
        qrsp4 = int16(reshape(qrsp4, [thsNtrl 1]));
        res.Ib4(repi,sampi) = info.calc_info_integer_c_int16_t(qrsp4, 16, qstm, 2, thsNtrl);
        res.Ib4mm(repi,sampi) = res.Ib4(repi,sampi) - 15./(2*thsNtrl*log(2));

        rqrsp8 = double(bin.eqpop(squeeze(thseeg(:,1)),8));
        dqrsp8 = double(bin.eqpop(squeeze(thseeg(:,2)),8));
        qrsp8 = numbase2dec([rqrsp8(:) dqrsp8(:)]',8);
        qrsp8 = int16(reshape(qrsp8, [thsNtrl 1]));
        res.Ib8(repi,sampi) = info.calc_info_integer_c_int16_t(qrsp8, 64, qstm, 2, thsNtrl);
        res.Ib8mm(repi,sampi) = res.Ib8(repi,sampi) - 63./(2*thsNtrl*log(2));
        
        crsp = copnorm(thseeg);
%         try
            res.Icop(repi,sampi) = info_gd(crsp, qstm, 2, true, false, false);
%         catch
%             res.Icop(repi,sampi) = NaN;
%         end
        if ~isfinite(res.Icop(repi,sampi))
            keyboard
        end
    end
end

%%
est = {'Ib2' 'Ib4' 'Ib8' 'Ib2mm' 'Ib4mm' 'Ib8mm' 'Icop'};
figure
hold all
for ei=1:length(est)
    plot(2.^log_samp_size, nanmean(res.(est{ei}),1));
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