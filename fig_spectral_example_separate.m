% spectral MI figure with 2 examples
% separate phase and power

%% Example 1

Nsamp = 500;
Ny = 2;
Nplt = 300;

%%
y = randi(Ny,[1 Nsamp]);

th = zeros(1,Nsamp);
r = zeros(1,Nsamp);

vmmu = [0 3.14];
exmu = 5*[1 1];
for yi=1:Ny
    idx = y==yi;
    thsN = sum(idx);
    th(idx) = vmrand(vmmu(yi),1,1,thsN);
%     r(idx) = exprnd(exmu(yi),1,thsN);
    r(idx) = chi2rnd(exmu(yi),1,thsN);
end
x = (r .* exp(1i.*th))';

figure
al = 0.2;
cols = {'r' 'b'};

subplot(4,4,1)
hold on
for yi=1:Ny
    idx = y==yi;
    thsx = x(idx);
    hs = scatter(real(thsx),imag(thsx),cols{yi},'filled');
    alpha(hs,al)
end
axis square

cx = [copnorm(real(x)) copnorm(imag(x))];
subplot(4,4,5)
hold on
mu = mean(cphs);
sig = cov(cphs);
xgg = linspace(-5,5,100);
ygg = linspace(-5,5,100);
[xg yg] = meshgrid(xgg,ygg);
z = mvnpdf([xg(:) yg(:)],mu,sig);
z = reshape(z,size(xg));
v = [0.01 0.01];
contour(xg,yg,z,v,'k')
for yi=1:Ny
    idx = y==yi;
    thsx = cx(idx,:);
    hs = scatter(thsx(:,1),thsx(:,2),cols{yi},'filled');
    alpha(hs,al)
    mu = mean(thsx);
    sig = cov(thsx);
    z = mvnpdf([xg(:) yg(:)],mu,sig);
    z = reshape(z,size(xg));
    contour(xg,yg,z,v,cols{yi})
end
axis square

pow = abs(x);
cpow = copnorm(pow);
subplot(4,4,2);
hold on
[f fxi] = ksdensity(pow,'support','positive');
plot(fxi,f,'k')
for yi=1:Ny
    idx = y==yi;
    ths = pow(idx);
    [f fxi] = ksdensity(ths,fxi,'support','positive');
    plot(fxi,f,cols{yi})
end

subplot(4,4,6);
hold on
[f fxi] = ksdensity(cpow);
plot(fxi,f,'k')
for yi=1:Ny
    idx = y==yi;
    ths = cpow(idx);
    [f fxi] = ksdensity(ths,fxi);
    plot(fxi,f,cols{yi})
end


xphs = x ./ abs(x);
subplot(4,4,3)
hold on
for yi=1:Ny
    idx = y==yi;
    ths = xphs(idx);
    hs = scatter(real(ths),imag(ths),cols{yi},'filled');
    alpha(hs,al);
end
axis square

cphs = copnorm([real(xphs) imag(xphs)]);
subplot(4,4,7)
hold on
for yi=1:Ny
    idx = y==yi;
    ths = cphs(idx,:);
    hs = scatter(ths(:,1),ths(:,2),cols{yi},'filled');
    alpha(hs,al);
end

mu = mean(cphs);
sig = cov(cphs);
xgg = linspace(-5,5,100);
ygg = linspace(-5,5,100);
[xg yg] = meshgrid(xgg,ygg);
z = mvnpdf([xg(:) yg(:)],mu,sig);
z = reshape(z,size(xg));
v = [0.01 0.01];
contour(xg,yg,z,v,'k')
% figure
% subplot(1,3,1);
% imagesc(xgg,ygg,z);colorbar
for yi=1:Ny
    idx = y==yi;
    ths = cphs(idx,:);
%     hs = scatter(ths(:,1),ths(:,2),cols{yi},'filled');
%     alpha(hs,al);
    mu = mean(ths);
    sig = cov(ths);
    z = mvnpdf([xg(:) yg(:)],mu,sig);
    z = reshape(z,size(xg));
    contour(xg,yg,z,v,cols{yi})
%     subplot(1,3,yi+1)
%     imagesc(xgg,ygg,z);colorbar
end
axis square


subplot(2,4,4)
Ispec = info_gd(cx, y-1, 2, true, false, false);
Ipow = info_gd(cpow, y-1, 2, true, false, false);
Iphs = info_gd(cphs, y-1, 2, true, false, false);
% Iphspow = info_gd([cphs cpow], y-1, 2, true, false,false)
bar([Ispec Ipow Iphs])
ylim([0 0.35])

% %%
y = randi(Ny,[1 Nsamp]);

th = zeros(1,Nsamp);
r = zeros(1,Nsamp);

vmmu = [0 0];
exmu = 5*[1 2];
for yi=1:Ny
    idx = y==yi;
    thsN = sum(idx);
    th(idx) = vmrand(vmmu(yi),0.3,1,thsN);
%     r(idx) = exprnd(exmu(yi),1,thsN);
    r(idx) = chi2rnd(exmu(yi),1,thsN);
end
x = (r .* exp(1i.*th))';

% figure
al = 0.2;
cols = {'r' 'b'};

subplot(4,4,9)
hold on
for yi=1:Ny
    idx = y==yi;
    thsx = x(idx);
    hs = scatter(real(thsx),imag(thsx),cols{yi},'filled');
    alpha(hs,al)
end
axis square

cx = [copnorm(real(x)) copnorm(imag(x))];
subplot(4,4,13)
hold on
mu = mean(cphs);
sig = cov(cphs);
xgg = linspace(-5,5,100);
ygg = linspace(-5,5,100);
[xg yg] = meshgrid(xgg,ygg);
z = mvnpdf([xg(:) yg(:)],mu,sig);
z = reshape(z,size(xg));
v = [0.01 0.01];
contour(xg,yg,z,v,'k')
for yi=1:Ny
    idx = y==yi;
    thsx = cx(idx,:);
    hs = scatter(thsx(:,1),thsx(:,2),cols{yi},'filled');
    alpha(hs,al)
    mu = mean(thsx);
    sig = cov(thsx);
    z = mvnpdf([xg(:) yg(:)],mu,sig);
    z = reshape(z,size(xg));
    contour(xg,yg,z,v,cols{yi})
end
axis square

pow = abs(x);
cpow = copnorm(pow);
subplot(4,4,10);
hold on
[f fxi] = ksdensity(pow,'support','positive');
plot(fxi,f,'k')
for yi=1:Ny
    idx = y==yi;
    ths = pow(idx);
    [f fxi] = ksdensity(ths,fxi,'support','positive');
    plot(fxi,f,cols{yi})
end

subplot(4,4,14);
hold on
[f fxi] = ksdensity(cpow);
plot(fxi,f,'k')
for yi=1:Ny
    idx = y==yi;
    ths = cpow(idx);
    [f fxi] = ksdensity(ths,fxi);
    plot(fxi,f,cols{yi})
end


xphs = x ./ abs(x);
subplot(4,4,11)
hold on
for yi=1:Ny
    idx = y==yi;
    ths = xphs(idx);
    hs = scatter(real(ths),imag(ths),cols{yi},'filled');
    alpha(hs,al);
end
axis square

cphs = copnorm([real(xphs) imag(xphs)]);
subplot(4,4,15)
hold on
for yi=1:Ny
    idx = y==yi;
    ths = cphs(idx,:);
    hs = scatter(ths(:,1),ths(:,2),cols{yi},'filled');
    alpha(hs,al);
end
axis square

mu = mean(cphs);
sig = cov(cphs);
xgg = linspace(-5,5,100);
ygg = linspace(-5,5,100);
[xg yg] = meshgrid(xgg,ygg);
z = mvnpdf([xg(:) yg(:)],mu,sig);
z = reshape(z,size(xg));
v = [0.01 0.01];
contour(xg,yg,z,v,'k')
% figure
% subplot(1,3,1);
% imagesc(xgg,ygg,z);colorbar
for yi=1:Ny
    idx = y==yi;
    ths = cphs(idx,:);
%     hs = scatter(ths(:,1),ths(:,2),cols{yi},'filled');
%     alpha(hs,al);
    mu = mean(ths);
    sig = cov(ths);
    z = mvnpdf([xg(:) yg(:)],mu,sig);
    z = reshape(z,size(xg));
    contour(xg,yg,z,v,cols{yi})
%     subplot(1,3,yi+1)
%     imagesc(xgg,ygg,z);colorbar
end

subplot(2,4,8)
Ispec = info_gd(cx, y-1, 2, true, false, false);
Ipow = info_gd(cpow, y-1, 2, true, false, false);
Iphs = info_gd(cphs, y-1, 2, true, false, false);
% Iphspow = info_gd([cphs cpow], y-1, 2, true, false,false)
bar([Ispec Ipow Iphs])
ylim([0 0.35])
