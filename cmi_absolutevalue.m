Ntrl = 100000;

x = 8*(rand(Ntrl,1)-0.5);
ax = abs(x);

y = randn(Ntrl,1) + x;
z = randn(Ntrl,1) + ax;

%%
Npt = 500;
idx = 1:Npt;
figure
subplot(321)
scatter(x(idx),ax(idx),'k','filled')
xlabel('x')
ylabel('ax')
subplot(323)
scatter(x(idx),y(idx),'k','filled')
xlabel('x')
ylabel('y')
subplot(324)
scatter(ax(idx),y(idx),'k','filled')
xlabel('ax')
ylabel('y')
subplot(325)
scatter(x(idx),z(idx),'k','filled')
xlabel('x')
ylabel('z')
subplot(326)
scatter(ax(idx),z(idx),'k','filled')
xlabel('ax')
ylabel('z')


%%
cx = copnorm(x);
cax = copnorm(ax);
cy = copnorm(y);
cz = copnorm(z);

% Iyx = info_gg(cy,cx,true,false,false);
% Iyax = info_gg(cy,cax,true,false,false);
% Izx = info_gg(cz,cx,true,false,false);
% Izax = info_gg(cz,cax,true,false,false);

Iyx = info_gg(cy,cx,true,false,false)
Iyax = info_gg(cy,cax,true,false,false)
Izx = info_gg(cz,cx,true,false,false)
Izax = info_gg(cz,cax,true,false,false)

IyxCax = cmi_ggg(cy,cx,cax,true,false,false)
IyaxCx = cmi_ggg(cy,cax,cx,true,false,false)
IzxCax = cmi_ggg(cz,cx,cax,true,false,false)
IzaxCx = cmi_ggg(cz,cax,cx,true,false,false)

%%
% Nbin = 3;
% bx = bin.eqpop(x,Nbin);
% bax = bin.eqpop(ax,Nbin);
% by = bin.eqpop(y,Nbin);
% bz = bin.eqpop(z,Nbin);
% 
% Iyx = info.calc_info(by,Nbin,bx,Nbin,Ntrl)
% Iyax = info.calc_info(by,Nbin,bax,Nbin,Ntrl)
% Izx = info.calc_info(bz,Nbin,bx,Nbin,Ntrl)
% Izax = info.calc_info(bz,Nbin,bax,Nbin,Ntrl)
% 
% IyxCax = info.calc_cmi(by,Nbin,bx,Nbin,bax,Nbin,Ntrl)
% IyaxCx = info.calc_cmi(by,Nbin,bax,Nbin,bx,Nbin,Ntrl)
% IzxCax = info.calc_cmi(bz,Nbin,bx,Nbin,bax,Nbin,Ntrl)
% IzaxCx = info.calc_cmi(bz,Nbin,bax,Nbin,bx,Nbin,Ntrl)

%%
Cyx = corr(y,x)
Cyax = corr(y,ax)
Czx = corr(z,x)
Czax = corr(z,ax)

%%

%%
[coef score] = princomp([ax y]);
figure
subplot(1,2 ,1)
hs = scatter(ax,y,'filled');
alpha(hs,0.2)
subplot(1,2,2)
hs = scatter(score(:,1),score(:,2),'filled');
alpha(hs,0.2)

