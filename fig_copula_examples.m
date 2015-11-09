
%% uncorrelated Gaussian copula example

Npoint = 1000;
x = mvnrnd([0 0],[1 0; 0 1],Npoint);
figure
subplot(1,2,1)
hs = scatter(x(:,1),x(:,2),16,'k','filled');
set(gca,'ticklength',2*get(gca,'ticklength'))
alpha(hs,al);
box on
axis square
l = 4;
xlim([-l l])
ylim([-l l])

subplot(1,2,2)
hs = scatter(ctransform(x(:,1)), ctransform(x(:,2)), 'k','filled');
set(gca,'ticklength',2*get(gca,'ticklength'))
alpha(hs,al)
axis square
box on
%% correlated Gaussian copula example
Npoint = 1000;
r = 0.8;
al = 0.2;

x = mvnrnd([0 0],[1 r; r 1],Npoint);
figure
subplot(1,2,1)
hs = scatter(x(:,1),x(:,2),16,'k','filled');
set(gca,'ticklength',2*get(gca,'ticklength'))
alpha(hs,al);
box on
axis square
l = 4;
xlim([-l l])
ylim([-l l])

subplot(1,2,2)
hs = scatter(ctransform(x(:,1)), ctransform(x(:,2)), 'k','filled');
set(gca,'ticklength',2*get(gca,'ticklength'))
alpha(hs,al)
axis square
box on

%% Gaussian-exponentialexample
Npoint = 1000;
mu = 1.5;
y = exprnd(mu, [Npoint 1]);
x = randn(Npoint,1) + y;

figure
subplot(1,3,1)
% scatter(y,x,'k','filled')
scatterhist(y,x,'Kernel','on','Color','k','Marker','.')
set(gca,'ticklength',2*get(gca,'ticklength'))
% % 
% axis square
% l = 4;
% xlim([0 14])
% ylim([-4 14])

%%
cx = ctransform(x);
cy = ctransform(y);
figure
subplot(1,2,1)
hs = scatter(cy, cx, 'k','filled');
set(gca,'ticklength',2*get(gca,'ticklength'))
alpha(hs,al)
axis square
box on
%%
subplot(1,3,3)
scatter(norminv(cy),norminv(cx),'k','filled')
axis square
%%
scatterhist(norminv(cy),norminv(cx),'Kernel','on','Color','k','Marker','.')
set(gca,'ticklength',2*get(gca,'ticklength'))


%% Discrete example - RAW
Npoint = 500;
al = 0.1;
y = randi(2,[Npoint 1]);
sh = 3;
x = (randn(Npoint,1)-(sh/2)) + sh*(y-1);
% x = copnorm(x)
Nscat = 500;
xl = 5;

lw = 2;
al = 0.1;
figure
subplot(3,1,2)
% [f fxi] = ksdensity(x);
fxi = linspace(-xl,xl,100);
% f = (normpdf(fxi)+normpdf(fxi-sh)+normpdf(fxi+sh))./3;
f = (normpdf(fxi+(sh/2))+normpdf(fxi-(sh/2)))./2;
full.f = f;
full.fxi = fxi;
% plot(fxi,normpdf(fxi),'k')
plot(fxi,f,'k','LineWidth',lw)
hold on
hs = scatter(x(1:Nscat),0.3+(randn(Nscat,1)/100),'k','filled');
alpha(hs,al)
xlim([-xl xl])
box off
plot(fxi, normpdf(fxi, mean(x),std(x)),'k:')

subplot(3,1,1)
hold on
cols = {'r' 'b'};
offsets = [0.5 0.6];
grp = [];
for i=1:2
    thsx = x(y==i);
%     [f fxi] = ksdensity(thsx);
    f = normpdf(fxi-sh*(i-1)+(sh/2));
    grp(i).f = f;
    grp(i).fxi = fxi;
    plot(fxi, f, cols{i},'LineWidth',lw)
    thsNpt = min(Nscat,length(thsx));
    hs = scatter(thsx(1:thsNpt), offsets(i)+(randn(thsNpt,1)/100),cols{i},'filled');
    alpha(hs,al)
end
xlim([-xl xl])

subplot(3,1,3)
hold on
plot(full.fxi,full.f,'k','LineWidth',lw);
fxi = linspace(-5,5,100);
% grpf = zeros(100,2);
% for i=1:2
%     thsx = x(y==i);
%     [f ~] = ksdensity(thsx-mean(thsx),fxi);
% %     f = normpdf(fxi);
%     plot(fxi,f,cols{i})
%     grpf(:,i) = f;
% end
% plot(fxi,mean(grpf,2),'--k','LineWidth',2)

xlim([-xl xl])
plot(fxi,normpdf(fxi),'--k','LineWidth',2)
plot(fxi, normpdf(fxi, mean(x),std(x)),'k:')

%% Discrete example - COPULA
Npoint = 500;
% y = randi(2,[Npoint 1]);
sh = 3;
% x = (randn(Npoint,1)-(sh/2)) + sh*(y-1);
x = copnorm(x);
Nscat = 500;
xl = 3.5;

lw = 2;
al = 0.1;
figure
subplot(3,1,2)
% [f fxi] = ksdensity(x);
fxi = linspace(-xl,xl,100);
% f = (normpdf(fxi)+normpdf(fxi-sh)+normpdf(fxi+sh))./3;
% f = (normpdf(fxi+(sh/2))+normpdf(fxi-(sh/2)))./2;
f = normpdf(fxi);
full.f = f;
full.fxi = fxi;
% plot(fxi,normpdf(fxi),'k')
plot(fxi,f,'k','LineWidth',lw)
hold on
hs = scatter(x(1:Nscat),0.5+(randn(Nscat,1)/75),'k','filled');
alpha(hs,al)
xlim([-xl xl])
box off
plot(fxi, normpdf(fxi, mean(x),std(x)),'k:')

subplot(3,1,1)
hold on
cols = {'r' 'b'};
offsets = [0.8 0.9];
grp = [];
for i=1:2
    thsx = x(y==i);
    [f ~] = ksdensity(thsx,fxi);
%     f = normpdf(fxi-sh*(i-1)+(sh/2));
    grp(i).f = f;
    grp(i).fxi = fxi;
    plot(fxi, f, cols{i},'LineWidth',lw)
    thsNpt = min(Nscat,length(thsx));
    hs = scatter(thsx(1:thsNpt), offsets(i)+(randn(thsNpt,1)/100),cols{i},'filled');
    alpha(hs,al)
end
xlim([-xl xl])

subplot(3,1,3)
hold on
plot(full.fxi,full.f,'k','LineWidth',lw);
fxi = linspace(-5,5,100);
grpf = zeros(100,2);
for i=1:2
    thsx = x(y==i);
    [f ~] = ksdensity(thsx-mean(thsx),fxi);
%     f = normpdf(fxi);
    plot(fxi,f,cols{i})
    grpf(:,i) = f;
end
% plot(fxi,mean(grpf,2),'--k','LineWidth',2)

xlim([-xl xl])
% plot(fxi,normpdf(fxi),'--k','LineWidth',2)
% plot(fxi, normpdf(fxi, mean(x),std(x)),'k:')
    
%% Discrete example - no effect
Npoint = 500;
y = randi(2,[Npoint 1]);
sh = 3;
% x = (randn(Npoint,1)-(sh/2)) + sh*(y-1);
x = randn(Npoint,1);
% x = copnorm(x)
Nscat = 500;
xl = 5;

lw = 2;
al = 0.2;
figure
subplot(3,1,2)
% [f fxi] = ksdensity(x);
fxi = linspace(-xl,xl,100);
% f = (normpdf(fxi)+normpdf(fxi-sh)+normpdf(fxi+sh))./3;
% f = (normpdf(fxi+(sh/2))+normpdf(fxi-(sh/2)))./2;
f = normpdf(fxi);
full.f = f;
full.fxi = fxi;
% plot(fxi,normpdf(fxi),'k')
plot(fxi,f,'k','LineWidth',lw)
hold on
hs = scatter(x(1:Nscat),0.3+(randn(Nscat,1)/100),'k','filled');
alpha(hs,al)
xlim([-xl xl])
box off
plot(fxi, normpdf(fxi, mean(x),std(x)),'k:')

subplot(3,1,1)
hold on
cols = {'r' 'b'};
offsets = [0.5 0.6];
grp = [];
for i=1:2
    thsx = x(y==i);
%     [f fxi] = ksdensity(thsx);
%     f = normpdf(fxi-sh*(i-1)+(sh/2));
    f = normpdf(fxi);
    grp(i).f = f;
    grp(i).fxi = fxi;
    plot(fxi, f, cols{i},'LineWidth',lw)
    thsNpt = min(Nscat,length(thsx));
    hs = scatter(thsx(1:thsNpt), offsets(i)+(randn(thsNpt,1)/100),cols{i},'filled');
    alpha(hs,al)
end
xlim([-xl xl])

subplot(3,1,3)
hold on
plot(full.fxi,full.f,'k','LineWidth',lw);
fxi = linspace(-5,5,100);
% grpf = zeros(100,2);
% for i=1:2
%     thsx = x(y==i);
%     [f ~] = ksdensity(thsx-mean(thsx),fxi);
% %     f = normpdf(fxi);
%     plot(fxi,f,cols{i})
%     grpf(:,i) = f;
% end
% plot(fxi,mean(grpf,2),'--k','LineWidth',2)

xlim([-xl xl])
plot(fxi,normpdf(fxi),'--k','LineWidth',2)
plot(fxi, normpdf(fxi, mean(x),std(x)),'k:')
    
% set(get(gcf,'Children'),'ticklength',2*get(gca,'ticklength'))



