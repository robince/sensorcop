

%% features
mu = [0 0];
C = 0.7;
sigma = [1 C; C 1];
x = linspace(-5,5,100);
y = linspace(-5,5,100);

[X Y] = meshgrid(x,y);
Z = mvnpdf([X(:) Y(:)],mu,sigma);
Z = reshape(Z,size(X));



%%
figure

imagesc(x,y,Z)
set(gca,'YDir','normal')
axis square



%%
mu = [0 0];
C = 0.6;
sigma = [1 C; C 1];
Npt = 50000;
% generate features
stim = mvnrnd(mu,sigma,Npt);
respA = randn(Npt,1) + stim(:,1);
respB = randn(Npt,1) + (stim(:,1)+stim(:,2))./2;

%%
figure
scatter(stim(:,1),stim(:,2),'k','filled')
xlabel('Stim A')
ylabel('Stim B')
title('Resp = Stim A + Noise')

%%
figure
xl = [-5 5];
Nplt = 500;
idx = randperm(Npt, Nplt);
subplot(2,2,1)
scatter(respA(idx), stim(idx,1), 'k','filled')
xlabel('Resp A')
ylabel('Stim A')
mi = info_gg(stim(:,1),respA, true, false, false);
cmi = cmi_ggg(stim(:,1),respA,stim(:,2), true, false, false);
title(sprintf('MI: %.3f  CMI: %.3f',mi,cmi));
xlim(xl)
ylim(xl)

subplot(2,2,3)
scatter(respA(idx), stim(idx,2), 'k','filled')
xlabel('Resp A')
ylabel('Stim B')
mi = info_gg(stim(:,2),respA, true, false, false);
cmi = cmi_ggg(stim(:,2),respA,stim(:,1), true, false, false);
title(sprintf('MI: %.3f  CMI: %.3f',mi,cmi));
xlim(xl)
ylim(xl)

subplot(2,2,2)
scatter(respB(idx), stim(idx,1), 'k','filled')
xlabel('Resp B')
ylabel('Stim A')
mi = info_gg(stim(:,1),respB, true, false, false);
cmi = cmi_ggg(stim(:,1),respB,stim(:,2), true, false, false);
title(sprintf('MI: %.3f  CMI: %.3f',mi,cmi));
xlim(xl)
ylim(xl)

subplot(2,2,4)
scatter(respB(idx), stim(idx,2), 'k','filled')
xlabel('Resp B')
ylabel('Stim B')
mi = info_gg(stim(:,2),respB, true, false, false);
cmi = cmi_ggg(stim(:,2),respB,stim(:,1), true, false, false);
title(sprintf('MI: %.3f  CMI: %.3f',mi,cmi));
xlim(xl)
ylim(xl)

%%
figure
subplot(1,2,1)
% scatter3(resp, stim(:,1), stim(:,2), 'k', 'filled')
plotmatrix(respB, [stim(:,1), stim(:,2)])
subplot(1,2,2)
plotmatrix(respA, [stim(:,1), stim(:,2)])