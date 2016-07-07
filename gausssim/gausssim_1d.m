% Gaussian simulation to look at MI value bias / MSE
% univariate gauss

cor = [0.2 0.4 0.6 0.8];
Ncor = length(cor);

log_samp_size = 3:10;
Nsamp = length(log_samp_size);

k = [2 3 4];
Nk = length(k);

Nrep = 500;

Icop = zeros(Nrep,Ncor,Nsamp);
Ib2 = zeros(Nrep,Ncor,Nsamp);
Ib4 = zeros(Nrep,Ncor,Nsamp);
Ib8 = zeros(Nrep,Ncor,Nsamp);
Ik1 = zeros(Nrep,Ncor,Nsamp,Nk);
Ik2 = zeros(Nrep,Ncor,Nsamp,Nk);


for si=1:Nsamp
    samps = 2.^log_samp_size(si)
    for ci=1:Ncor
        c = cor(ci);
        
        parfor ri=1:Nrep
            
            % generate data
            d = mvnrnd([ 0 0], [1 c; c 1], samps);
            x = d(:,1);
            y = d(:,2);
            
            Icop(ri,ci,si) = gcmi_cc(x,y);
            
            Nbin = 2;
            bx = int16(bin.eqpop(x,Nbin));
            by = int16(bin.eqpop(y,Nbin));
            Ib2(ri,ci,si) = info.calc_info_integer_c_int16_t(bx,Nbin,by,Nbin,samps);
            Nbin = 4;
            bx = int16(bin.eqpop(x,Nbin));
            by = int16(bin.eqpop(y,Nbin));
            Ib4(ri,ci,si) = info.calc_info_integer_c_int16_t(bx,Nbin,by,Nbin,samps);
            Nbin = 8;
            bx = int16(bin.eqpop(x,Nbin));
            by = int16(bin.eqpop(y,Nbin));
            Ib8(ri,ci,si) = info.calc_info_integer_c_int16_t(bx,Nbin,by,Nbin,samps);
            
            for ki=1:Nk
                [I1 I2] = KraskovMI(x,y,k(ki));
                Ik1(ri,ci,si,ki) = I1;
                Ik2(ri,ci,si,ki) = I2;
            end
        end
    end
end

Ik1 = Ik1 ./ log(2);
Ik2 = Ik2 ./ log(2);

%%
save('1dres_500.mat')
%%
% c = 0.2;
% Itrue = (log(2*pi*exp(1)) - 1 - log(2*pi) - 0.5*log(det([ 1 c; c 1]))) ./ log(2)
binbias = 1./((2.^log_samp_size)*2*log(2));
Itrue = zeros(1,Ncor);
for ci=1:Ncor
    c = cor(ci);
    Itrue(ci) = (log(2*pi*exp(1)) - 1 - log(2*pi) - 0.5*log(det([ 1 c; c 1])))./ log(2);
end
%%

figure
subplot(2,4,1)
hold all
sidx = 2:Nsamp;
plot(log_samp_size(sidx), squeeze(mean(squeeze(Icop(:,1,sidx)),1))')
xlim([3.5 10.5])
plot(log_samp_size(sidx), squeeze(mean(squeeze(Ib2(:,1,sidx)),1))' - binbias(sidx)')
plot(log_samp_size(sidx), squeeze(mean(squeeze(Ib4(:,1,sidx)),1))' - 9*binbias(sidx)')
plot(log_samp_size(sidx), squeeze(mean(squeeze(Ik1(:,1,sidx,1)),1))')
plot(log_samp_size(sidx), squeeze(mean(squeeze(Ik1(:,2,sidx,1)),1))')
% plot(log_samp_size(sidx), squeeze(mean(squeeze(Ik1(:,3,sidx,1)),1))')
ylim([0 0.25])

%% MEAN / IQR
sidx = 2:Nsamp;
x = log_samp_size(sidx);

figure

for ci=1:4

subplot(2,4,ci)
pd = squeeze(Icop(:,ci,sidx));
errorbar(x-0.1, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));
hold all
% 
pd =bsxfun(@minus,squeeze(Ib2(:,ci,sidx)),binbias(sidx));
errorbar(x-0.05, median(pd), median(pd)-prctile(pd,25),prctile(pd,75)-median(pd));

pd = bsxfun(@minus,squeeze(Ib4(:,ci,sidx)),9*binbias(sidx));
errorbar(x+0.05, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));

% pd = bsxfun(@minus,squeeze(Ib8(:,ci,sidx)),49*binbias(sidx));
% errorbar(x, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));

% pd = squeeze(Ik1(:,ci,sidx,1));
% errorbar(x+0.1, median(pd), median(pd)-prctile(pd,25),prctile(pd,75)-median(pd));
% 
pd = squeeze(Ik1(:,ci,sidx,2));
errorbar(x+0.1, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));

% pd = squeeze(Ik1(:,ci,sidx,3));
% errorbar(x+0.3, median(pd), median(pd)-prctile(pd,25),prctile(pd,75)-median(pd));

% subplot(2,4,2)
% pd = squeeze(Icop(:,ci,sidx));
% errorbar(x, mean(pd), var(pd)/sqrt(size(pd,1))  );

xlim([3.5 10.5])
ylim([-0.15 1])
hline(Itrue(ci),'k:')
title(sprintf('c = %.1f', cor(ci)))
end
subplot(2,4,1)
legend('GCMI','2 bin','4 bin','kNN')
ylabel('MI (bits)')

% %% mean square error
for ci=1:4

subplot(2,4,4+ci)
pd = squeeze(Icop(:,ci,sidx));
pd = (pd - Itrue(ci)).^2;
sem = std(pd) ./ sqrt(size(pd,1));
errorbar(x-0.1, mean(pd), sem./2);
hold all

pd =bsxfun(@minus,squeeze(Ib2(:,ci,sidx)),binbias(sidx));
pd = (pd - Itrue(ci)).^2;
sem = std(pd) ./ sqrt(size(pd,1));
errorbar(x-0.05, mean(pd), sem./2);

pd = bsxfun(@minus,squeeze(Ib4(:,ci,sidx)),9*binbias(sidx));
pd = (pd - Itrue(ci)).^2;
sem = std(pd) ./ sqrt(size(pd,1));
errorbar(x+0.05, mean(pd), sem./2);
%
pd = squeeze(Ik1(:,ci,sidx,2));
pd = (pd - Itrue(ci)).^2;
sem = std(pd) ./ sqrt(size(pd,1));
errorbar(x+0.1, mean(pd), sem./2);

xlim([3.5 10.5])
xlabel('log_2 samples')
% ylim([-0.15 1])
% hline(Itrue(ci),'k:')
end
subplot(2,4,1)
ylabel('MSE')


%% MEAN / IQR EFFECT OF K
sidx = 2:Nsamp;
x = log_samp_size(sidx);

figure

for ci=1:4

subplot(2,4,ci)
% pd = squeeze(Icop(:,ci,sidx));
% errorbar(x-0.1, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));
hold all
% 
% pd =bsxfun(@minus,squeeze(Ib2(:,ci,sidx)),binbias(sidx));
% errorbar(x-0.05, median(pd), median(pd)-prctile(pd,25),prctile(pd,75)-median(pd));

% pd = bsxfun(@minus,squeeze(Ib4(:,ci,sidx)),9*binbias(sidx));
% errorbar(x+0.05, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));

% pd = bsxfun(@minus,squeeze(Ib8(:,ci,sidx)),49*binbias(sidx));
% errorbar(x, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));

% pd = squeeze(Ik1(:,ci,sidx,1));
% errorbar(x+0.1, median(pd), median(pd)-prctile(pd,25),prctile(pd,75)-median(pd));
% 
pd = squeeze(Ik1(:,ci,sidx,1));
errorbar(x+0.1, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));
pd = squeeze(Ik1(:,ci,sidx,2));
errorbar(x+0.1, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));
pd = squeeze(Ik1(:,ci,sidx,3));
errorbar(x+0.1, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));
% pd = squeeze(Ik1(:,ci,sidx,3));
% errorbar(x+0.3, median(pd), median(pd)-prctile(pd,25),prctile(pd,75)-median(pd));

% subplot(2,4,2)
% pd = squeeze(Icop(:,ci,sidx));
% errorbar(x, mean(pd), var(pd)/sqrt(size(pd,1))  );

xlim([3.5 10.5])
ylim([-0.15 1])
hline(Itrue(ci),'k:')
end
subplot(2,4,1)
legend('2','3','4')
%%
figure
boxplot(squeeze(Icop(:,4,:)),'PlotStyle','compact','Symbol','','Whisker',0)
ylim([0 1])
hold all
% boxplot(squeeze(Ik1(:,4,:,3)),'PlotStyle','compact','Symbol','','Whisker',0)
%%

c = 0.2;
d = mvnrnd([ 0 0], [1 c; c 1], 5000);
x = d(:,1);
y = abs(d(:,2));
I = gcmi_cc(x,y)
[I1 I2] = KraskovMI(x,y,3)


%%
figure
% scatter(ctransform(x),ctransform(y),'filled')
% scatter(copnorm(x),copnorm(y),'filled')
