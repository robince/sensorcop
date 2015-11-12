
Nsamp = 10000;

Ny = 2;
y = randi(Ny,[1 Nsamp]);

th = zeros(1,Nsamp);
r = zeros(1,Nsamp);


vmmu = [0 3.14];
% vmmu = [0 2];
exmu = 10*[1.5 2];
for yi=1:Ny
    idx = y==yi;
    thsN = sum(idx);
    th(idx) = vmrand(vmmu(yi),1,1,thsN);
%     r(idx) = exprnd(exmu(yi),1,thsN);
    r(idx) = chi2rnd(exmu(yi),1,thsN);
end
x = r .* exp(i.*th);

%%
figure
subplot(2,4,1)
hold on
cols = {'r' 'b'};
for yi=1:Ny
    idx = y==yi;
    thsx = x(idx);
    scatter(real(thsx),imag(thsx),cols{yi},'filled')
end
% 
% subplot(2,4,5)
% hold on
% cols = {'r' 'b'};
% for yi=1:Ny
%     idx = y==yi;
%     thsx = x(idx);
%     scatter(ctransform(real(thsx)),ctransform(imag(thsx)),cols{yi},'filled')
% end


subplot(2,4,5)
hold on
cx = copnorm([real(x); imag(x)]');
for yi=1:Ny
    idx = y==yi;
    thsx = cx(idx,:);
    scatter(thsx(:,1),thsx(:,2),cols{yi},'filled')
end

Ifull = info_gd(cx,y-1,2,true,false,false)

subplot(2,4,2)
hold on
cpow = (abs(x)');
% fxi = linspace(0,3,100);
[f fxi] = ksdensity(cpow,'support','positive');
plot(fxi,f,'k')
for yi=1:Ny
    idx = y==yi;
    ths = cpow(idx);
    [f fxi] = ksdensity(ths,fxi,'support','positive');
    plot(fxi,f,cols{yi})
end
% xlim([0 3])


subplot(2,4,6)
hold on
cpow = copnorm(abs(x)');
[f fxi] = ksdensity(cpow);
plot(fxi,f,'k')
for yi=1:Ny
    idx = y==yi;
    ths = cpow(idx);
    [f fxi] = ksdensity(ths);
    plot(fxi,f,cols{yi})
end


xphs = x ./ abs(x);
subplot(2,4,3)
hold on
for yi=1:Ny
    idx = y==yi;
    ths = xphs(idx);
    scatter(real(ths),imag(ths),cols{yi},'filled')
end

cphs = copnorm([real(xphs); imag(xphs)]');
subplot(2,4,7)
hold on
for yi=1:Ny
    idx = y==yi;
    ths = cphs(idx,:);
    scatter(ths(:,1),ths(:,2),cols{yi},'filled')
end

%%
cphs = copnorm([real(exp(i*th)); imag(exp(i*th))]');
% rndr = rand(1,Nsamp);
% rndr = exprnd(1,[1 Nsamp]);
% cphs = copnorm([real(rndr.*exp(i*th)); imag(rndr.*exp(i*th))]');
Iphs = info_gd(cphs,y-1,2,true,false,false)

cpow = copnorm(r');
Ipow = info_gd(cpow,y-1,2,true,false,false)

% Ireal = info_gd(cx(:,1),y-1,2,true,false,false)
% Iimag = info_gd(cx(:,2),y-1,2,true,false,false)

Iphspow = info_gd([cphs cpow],y-1,2,true,false,false)
Isum = Iphs+Ipow

Itest = info_gd([cphs randn(size(cpow))],y-1,2,true,false,false)

% ISSUE WITH JOINT PHS+POW LESS THAN SUM
% SPECIFIC TO GD? (OR ALSO WITH PARAMETRIC GG)
%%
subplot(1,4,4)
bar([Ifull Iphs Ipow Iphspow])

%%
figure
scatter3(cphs(:,1),cphs(:,2),cpow,'filled')

%%
% TWO EXAMPLES - one phase modulated / one power modulated