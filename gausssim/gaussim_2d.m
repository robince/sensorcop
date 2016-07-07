% Gaussian simulation to look at MI value bias / MSE
% bivariate gauss

cor = [0.2 0.4 0.6];
Ncor = length(cor);

log_samp_size = 3:10;
Nsamp = length(log_samp_size);

k = 3;
Nk = length(k);

Nrep = 500;

Icop = zeros(Nrep,Ncor,Nsamp);
Ib2 = zeros(Nrep,Ncor,Nsamp);
Ib4 = zeros(Nrep,Ncor,Nsamp);
% Ib8 = zeros(Nrep,Ncor,Nsamp);
Ik1 = zeros(Nrep,Ncor,Nsamp,Nk);
Ik2 = zeros(Nrep,Ncor,Nsamp,Nk);


for si=1:Nsamp
    samps = 2.^log_samp_size(si)
    for ci=1:Ncor
        c = cor(ci);
        
        parfor ri=1:Nrep
            
            % generate data
            d = mvnrnd([0 0 0], [1 c c; c 1 0.1; c 0.1 1], samps);
            x = d(:,1);
            y = d(:,2:3);
            
            Icop(ri,ci,si) = gcmi_cc(x,y);
            
            Nbin = 2;
            bx = int16(bin.eqpop(x,Nbin));
            by1 = double(bin.eqpop(y(:,1),Nbin));
            by2 = double(bin.eqpop(y(:,2),Nbin));
            by = int16(numbase2dec([by1 by2]',2)');
            Ib2(ri,ci,si) = info.calc_info_integer_c_int16_t(bx,Nbin,by,Nbin.^2,samps);
            Nbin = 4;
            bx = int16(bin.eqpop(x,Nbin));
            by1 = double(bin.eqpop(y(:,1),Nbin));
            by2 = double(bin.eqpop(y(:,2),Nbin));
            by = int16(numbase2dec([by1 by2]',Nbin)');
            Ib4(ri,ci,si) = info.calc_info_integer_c_int16_t(bx,Nbin,by,Nbin.^2,samps);
%             Nbin = 8;
%             bx = int16(bin.eqpop(x,Nbin));
%             by1 = double(bin.eqpop(y(:,1),Nbin));
%             by2 = double(bin.eqpop(y(:,2),Nbin));
%             by = int16(numbase2dec([by1 by2]',2)');
%             Ib8(ri,ci,si) = info.calc_info_integer_c_int16_t(bx,Nbin,by,Nbin.^2,samps);
            
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
save('2dres_500.mat')
%%
binbias = 1./((2.^log_samp_size)*2*log(2));
Itrue = zeros(1,Ncor);
for ci=1:Ncor
    c = cor(ci);
    Hx = 0.5*log(2*pi*exp(1));
    Hy = 1+log(2*pi) + 0.5*log(det([1 0.1; 0.1 1]));
    Hxy = (3/2)*(1+log(2*pi)) + 0.5*log(det([1 c c; c 1 0.1; c 0.1 1]));
    Itrue(ci) = (Hx + Hy - Hxy)./ log(2);
end


%%

sidx = 2:Nsamp;
x = log_samp_size(sidx);

figure

for ci=1:3

subplot(2,3,ci)
pd = squeeze(Icop(:,ci,sidx));
errorbar(x-0.1, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));
hold all
% 
pd =bsxfun(@minus,squeeze(Ib2(:,ci,sidx)),3*binbias(sidx));
errorbar(x-0.05, median(pd), median(pd)-prctile(pd,25),prctile(pd,75)-median(pd));

pd = bsxfun(@minus,squeeze(Ib4(:,ci,sidx)),45*binbias(sidx));
errorbar(x+0.05, mean(pd), mean(pd)-prctile(pd,25),prctile(pd,75)-mean(pd));

% 
pd = squeeze(Ik1(:,ci,sidx,1));
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
subplot(2,3,1)
legend('GCMI','2 bin','4 bin','kNN')
ylabel('MI (bits)')


% %% mean square error
for ci=1:3

subplot(2,3,3+ci)
pd = squeeze(Icop(:,ci,sidx));
pd = (pd - Itrue(ci)).^2;
sem = std(pd) ./ sqrt(size(pd,1));
errorbar(x-0.1, mean(pd), sem./2);
hold all

pd =bsxfun(@minus,squeeze(Ib2(:,ci,sidx)),binbias(sidx));
pd = (pd - Itrue(ci)).^2;
sem = std(pd) ./ sqrt(size(pd,1));
errorbar(x-0.05, mean(pd), sem./2);

% pd = bsxfun(@minus,squeeze(Ib4(:,ci,sidx)),9*binbias(sidx));
% pd = (pd - Itrue(ci)).^2;
% sem = std(pd) ./ sqrt(size(pd,1));
% errorbar(x+0.05, mean(pd), sem./2);
plot(0,0)
%
pd = squeeze(Ik1(:,ci,sidx));
pd = (pd - Itrue(ci)).^2;
sem = std(pd) ./ sqrt(size(pd,1));
errorbar(x+0.1, mean(pd), sem./2);

xlim([3.5 10.5])
xlabel('log_2 samples')
% ylim([0 1])
% hline(Itrue(ci),'k:')
end
subplot(2,3,2)
ylabel('MSE')

