Nsamp = 200;
Ny = 2;

nodes = [0 pi/2 pi 3*pi/2]+pi/8;
yoffset = [0 pi/4];

%%
th = {};
y = {};
for yi=1:Ny
    for ni=1:length(nodes)
        thsth = vmrand(nodes(ni)+yoffset(yi), 15, 1, Nsamp);
        thsy = (yi-1)*ones(1,Nsamp);
        th{end+1} = thsth;
        y{end+1} = thsy;
    end
end

th = cell2mat(th);
y = cell2mat(y);
r = rand(size(th));
x = r .* exp(1i.*th);

%%
figure
subplot(121)
idx = y==0;
scatter(real(x(idx)),imag(x(idx)),'r','filled')
axis square
subplot(122)
idx = y==1;
scatter(real(x(idx)),imag(x(idx)),'b','filled')
axis square

%%
figure
al = 0.2;
cols = {'r' 'b'};

xphs = (x ./ abs(x)).';
subplot(2,1,1)
hold on
for yi=1:Ny
    idx = y==(yi-1);
    ths = xphs(idx);
    hs = scatter(real(ths),imag(ths),cols{yi},'filled');
    alpha(hs,al);
end
axis square

cphs = copnorm([real(xphs) imag(xphs)]);
subplot(2,1,2)
hold on
for yi=1:Ny
    idx = y==(yi-1);
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

for yi=1:Ny
    idx = y==(yi-1);
    ths = cphs(idx,:);
    mu = mean(ths);
    sig = cov(ths);
    z = mvnpdf([xg(:) yg(:)],mu,sig);
    z = reshape(z,size(xg));
    contour(xg,yg,z,v,cols{yi})
end
axis square

Hfull = ent_g(cphs,true,false)
Hclass = [ent_g(cphs(y==0,:),true,false) ent_g(cphs(y==1,:),true,false)]
Hcond = mean(Hclass)
I = Hfull - Hcond
