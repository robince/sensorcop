% continuous

figure
subplot(1,3,1)
x = linspace(-10,10,100);
sd = 1;
plot(x, normpdf(x, 0, sd),'k','LineWidth',2);
box off
H = 0.5 * log(2*pi*exp(1)*sd.^2) / log(2)

subplot(1,3,2)
sd = 4;
plot(x, normpdf(x, 0, sd),'k','LineWidth',2);
box off
H = 0.5 * log(2*pi*exp(1)*sd.^2) / log(2)


subplot(1,3,3)
plot(x, normpdf(x, -5, 1) + normpdf(x, 5, 1),'k','LineWidth',2);
box off
set(get(gcf,'Children'), 'YLim', [0 0.45])

%%
Nsamp = 400000;
x = [ randn(1,Nsamp)-5 randn(1,Nsamp)+5];
var(x)
%%
i = linspace(-15,15,1000);
p = ksdensity(x,i);
% -sum(p.*log2(p))
%%
-trapz(i,p.*log2(p))


%% discrete
figure
Nbin = 5;
subplot(1,3,1)
p = ones(Nbin,1);
p = p ./ sum(p);
bar((1:Nbin),p)
box off
xlim([0 Nbin+1])
ent(p)

subplot(1,3,2)
p = [0.5 3 6 3 0.5];
p = p ./ sum(p);
bar((1:Nbin),p)
box off
xlim([0 Nbin+1])
ent(p)

subplot(1,3,3)
p = [0.1 0.1 4 0.1  0.1];
p = p ./ sum(p);
bar((1:Nbin),p)
box off
xlim([0 Nbin+1])
ent(p)
set(get(gcf,'Children'), 'YLim', [0 1])


