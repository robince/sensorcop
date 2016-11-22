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

log_samp_size = 4:10;
Nsamp = length(log_samp_size);

Nrep = 500;
Ib2 = zeros(Nrep, Nsamp);
Ib4 = zeros(Nrep, Nsamp);
Ib2mm = zeros(Nrep, Nsamp);
Ib4mm = zeros(Nrep, Nsamp);
Icop = zeros(Nrep, Nsamp);
Ik = zeros(Nrep,Nsamp);

Nthread = 4;
for sampi=1:Nsamp
    thsNtrl = 2.^(log_samp_size(sampi))
    parfor repi=1:Nrep
        % jackknife (without replacement)
        idx = randperm(Ntrl, thsNtrl);
%         % bootstrap (with replacement)
%         idx = randi(Ntrl,thsNtrl,1);
        thsstim = stim(idx);
        ss = sum(thsstim==0);
        while (ss<4 || (thsNtrl-ss)<4)
%             idx = randperm(Ntrl, thsNtrl);
            idx = randi(Ntrl,thsNtrl,1);    
            thsstim = stim(idx);
            ss = sum(thsstim==0);
        end
        thseeg = rspdat(idx,:,ti,chi);
        qstm = int16(thsstim);
        % add jitter to avoid equal elements
        thseeg = thseeg + 1e-6*randn(size(thseeg));
        % true values
        rqrsp2 = double(bin.eqpop(squeeze(thseeg(:,1)),2));
        dqrsp2 = double(bin.eqpop(squeeze(thseeg(:,2)),2));
        qrsp2 = numbase2dec([rqrsp2(:) dqrsp2(:)]',2);
        qrsp2 = int16(reshape(qrsp2, [thsNtrl 1]));
        Ib2(repi,sampi) = info.calc_info_integer_c_int16_t(qrsp2, 4, qstm, 2, thsNtrl);
        Ib2mm(repi,sampi) = Ib2(repi,sampi) - 3./(2*thsNtrl*log(2));

        rqrsp4 = double(bin.eqpop(squeeze(thseeg(:,1)),4));
        dqrsp4 = double(bin.eqpop(squeeze(thseeg(:,2)),4));
        qrsp4 = numbase2dec([rqrsp4(:) dqrsp4(:)]',4);
        qrsp4 = int16(reshape(qrsp4, [thsNtrl 1]));
        Ib4(repi,sampi) = info.calc_info_integer_c_int16_t(qrsp4, 16, qstm, 2, thsNtrl);
        Ib4mm(repi,sampi) = Ib4(repi,sampi) - 15./(2*thsNtrl*log(2));

        crsp = copnorm(thseeg);
        Icop(repi,sampi) = info_gd(crsp, qstm, 2, true, false, false);
        
        % Kraskov
        k = 3;
        Hunc = KraskovEntropyV2(thseeg,k);
        w = zeros(1,2);
        Hcond = zeros(1,2);
        % class 0
        cidx = qstm==0;
        w(1) = sum(cidx) ./ thsNtrl;
        Hcond(1) = KraskovEntropyV2(thseeg(cidx,:),k);
        % class 1
        cidx = qstm==1;
        w(2) = sum(cidx) ./ thsNtrl;
        Hcond(2) = KraskovEntropyV2(thseeg(cidx,:),k);
        Ik(repi,sampi) = (Hunc - sum(w .* Hcond))/log(2);
    end
end
%%
clear dat rspdat drsp stim
save('2dbias_jack500.mat')
%%
est = { 'Ib2mm' 'Ib4mm'  'Icop' 'Ik'};
res = load('2dbias_boot500.mat');
figure
hold all
for ei=1:length(est)
    plot(log_samp_size, nanmean(res.(est{ei}),1));
end
legend(est)
%%
% est = {'Ib2' 'Ib4' 'Ib8' 'Icop'};
% figure
% hold all
% for ei=1:length(est)
%     plot(log_samp_size, nanstd(res.(est{ei}),1));
% end
% legend(est)
% %%
% est = {'Ib4'  'Icop'};
% figure
% hold all
% for ei=1:length(est)
%     errorbar(2.^log_samp_size, nanmean(res.(est{ei}),1),nanstd(res.(est{ei}),[],1));
% end
% legend(est)

%%
% load('2dbias_boot500.mat');
load('2dbias_jack500.mat');
sidx = 1:Nsamp;
x = log_samp_size(sidx);

figure

subplot(1,1,1)
pd = squeeze(Icop(:,sidx));
errorbar(x-0.1, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));
hold all

pd = squeeze(Ib2mm(:,sidx));
errorbar(x-0.033, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));

pd = squeeze(Ib4mm(:,sidx));
errorbar(x+0.033, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));

pd = squeeze(Ik(:,sidx));
errorbar(x+0.1, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));

xlim([3.5 10.5])
% ylim([-0.15 1])

% legend('GCMI','2 bin','4 bin','kNN')
ylabel('MI (bits)')
xlabel('log_2 samples')
% subplot(2,1,2)
