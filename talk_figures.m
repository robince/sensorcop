%% entropy examples
yl = [0 0.4];
Ns = 100000;
x = 10*randn(Ns,1);
[f,xi] = ksdensity(x);
figure;
area(xi,f)
ylim(yl)
axis off
box off

x = randn(Ns,1);
[f,xi] = ksdensity(x, xi);
figure
area(xi,f)
ylim(yl)
axis off
box off


x = [(randn(Ns/2,1) - 20); (randn(Ns/2,1)+20)];
[f,xi] = ksdensity(x, xi);
figure
area(xi,f)
ylim(yl)
axis off
box off


%% copula normalisation

x = (rand(Ns,1)*5)-2.5;
cx = copnorm(x);

figure
hold on
[f,xi] = ksdensity(cx,'function','cdf');
plot(xi,f,'m','LineWidth',2)
[f,xi] = ksdensity(x,'function', 'cdf');
plot(xi,f,'b','LineWidth',2)
box off
legend('Normal cdf', 'Empirical cdf','Location','SouthEast')
%%
figure
subplot(1,2,2)
[f,xi] = ksdensity(cx);
h = area(xi,f);
set(h,'FaceColor','m')
axis off
box off

subplot(1,2,1)
f = zeros(size(f));
f( xi>0 & xi<1 ) = 1;
h = area(xi,f);
set(h,'FaceColor','b')
axis off
box off