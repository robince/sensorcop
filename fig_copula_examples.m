
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
