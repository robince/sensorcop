load(fullfile(data_dir,'sub2_dat.mat'));
Fs = 50;

[Npln Nt] = size(plndat);

% pln sum
Nch = Npln / 2;
if round(Nch) ~= Nch
    error('Odd number of planar channels')
end

megH = plndat(1:Nch,:)';
megV = plndat((Nch+1):end,:)';
megdH = gradient_dim1(megH);
megdV = gradient_dim1(megV);
meg2d = permute(cat(3,megH,megV),[1 3 2]);
meg4d = permute(cat(3,megH,megV,megdH,megdV),[1 3 2]);

%%
chi = 231;
delay = 7;

log_samp_size = [8:13 14.45];
Nsamp = length(log_samp_size);

Nrep = 500;
Ib2 = zeros(Nrep, Nsamp);
Ib4 = zeros(Nrep, Nsamp);
Ib2mm = zeros(Nrep, Nsamp);
Ib4mm = zeros(Nrep, Nsamp);
Icop = zeros(Nrep, Nsamp);
Ik = zeros(Nrep,Nsamp);

%%
% align with delay
meg = meg2d((1+delay):end,:,chi);
spc = fltspc(1:(end-delay));
% normalise for scale (for adding jitter)
meg = bsxfun(@rdivide,meg,std(meg));
spc = spc ./ std(spc);
availsamp = length(spc);


Nthread = 4;
for sampi=Nsamp%1:(Nsamp-1)
    thsNtrl = round(2.^(log_samp_size(sampi)))
    parfor repi=1:Nrep
%         try
            %         idx = randperm(Ntrl, thsNtrl);
            
%             % bootstrap (with replacement)
%             idx = stationaryBB((1:availsamp)',thsNtrl,1,250);
            % jackknife (without replacement)
            idx = nonoverlappingBB((1:availsamp)',thsNtrl,1,250);
            thsspc = spc(idx);
            thsmeg = meg(idx,:);
            
            % add jitter to avoid equal elements
            thsspc = thsspc + 1e-6*randn(size(thsspc));
            thsmeg = thsmeg + 1e-6*randn(size(thsmeg));
            
            % true values
            rqrsp2 = double(bin.eqpop(squeeze(thsmeg(:,1)),2));
            dqrsp2 = double(bin.eqpop(squeeze(thsmeg(:,2)),2));
            qrsp2 = numbase2dec([rqrsp2(:) dqrsp2(:)]',2);
            qrsp2 = int16(reshape(qrsp2, [thsNtrl 1]));
            qstm = int16(bin.eqpop(thsspc,2));
            Ib2(repi,sampi) = info.calc_info_integer_c_int16_t(qrsp2, 4, qstm, 2, thsNtrl);
            Ib2mm(repi,sampi) = Ib2(repi,sampi)     - 3./(2*thsNtrl*log(2));
            
            rqrsp4 = double(bin.eqpop(squeeze(thsmeg(:,1)),4));
            dqrsp4 = double(bin.eqpop(squeeze(thsmeg(:,2)),4));
            qrsp4 = numbase2dec([rqrsp4(:) dqrsp4(:)]',4);
            qrsp4 = int16(reshape(qrsp4, [thsNtrl 1]));
            qstm = int16(bin.eqpop(thsspc,4));
            Ib4(repi,sampi) = info.calc_info_integer_c_int16_t(qrsp4, 16, qstm, 4, thsNtrl);
            Ib4mm(repi,sampi) = Ib4(repi,sampi) - 45./(2*thsNtrl*log(2));
            
            Icop(repi,sampi) = gcmi_cc(thsspc,thsmeg);
            
            % Kraskov
            k = 3;
            [I1] = KraskovMI(thsmeg,thsspc,k);
            Ik(repi,sampi) = I1 ./ log(2);
%         catch
%             disp(sprintf('Errored in rep %d',repi));
%             Ik(repi,sampi) = NaN;
%             Icop(repi,sampi) = NaN;
%             Ib2(repi,sampi) = NaN;
%             Ib2mm(repi,sampi) = NaN;
%             Ib4(repi,sampi) = NaN;
%             Ib4mm(repi,sampi) = NaN;
%         end
    end
end
%%
% save('megbias_boot500.mat','Ik','Icop','Ib2','Ib2mm','Ib4','Ib4mm','k','Nsamp','log_samp_size','chi','delay','Nrep')
save('megbias_jack500.mat','Ik','Icop','Ib2','Ib2mm','Ib4','Ib4mm','k','Nsamp','log_samp_size','chi','delay','Nrep')

%%
load('megbias_jack500.mat')


%%
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

% xlim([3.5 10.5])
% ylim([-0.15 1])

% legend('GCMI','2 bin','4 bin','kNN')
ylabel('MI (bits)')
xlabel('log_2 samples')
% subplot(2,1,2)