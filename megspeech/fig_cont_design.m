
lag = -34;
l = abs(lag);
x = fltspc(1:(end-l));
y = fltmeg((l+1):end);
% qspc = bin.eqpop(x,2);
qspc = zeros(size(x));
qspc(x>0.02) = 1;
Nplt = 5000;
idx = 50000+(1:Nplt);
time = (1:Nplt) ./ Fs;

%%
figure

thsx = x(idx);
thsqx = qspc(idx);

cut = max(x(qspc==0));
% 
% plttime = time;
% pltx = thsx;
% plttime(thsqx~=0) = NaN;
% pltx(thsqx~=0) = NaN;
% plot(plttime, pltx,'b')
% hold on
% plttime = time;
% pltx = thsx;
% plttime(thsqx~=1) = NaN;
% pltx(thsqx~=1) = NaN;
% plot(plttime, pltx,'r')
splitcolorplot(time,thsx',cut,'r','b')
box off
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02 0.025])
ylim([0 0.08])

pos = [680   827   528   148];
set(gcf,'Pos',pos);

%%
figure
plot(time,y(idx))
box off
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02 0.025])
pos = [680   827   528   148];
set(gcf,'Pos',pos);
%%
figure
kxi = linspace(-1e-12, 1e-12, 100);
[k0,~] = ksdensity(y(qspc==0),kxi);
[k1,~] = ksdensity(y(qspc==1),kxi);

figure
plot(kxi,k0, 'b','LineWidth',2);
hold on
plot(kxi,k1, 'r','LineWidth',2);
box off
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02 0.025])


%% random allignment

Npt = 1000;
Fs = 100;
x = randn(Npt,1);
y = randn(Npt,1);

figure
subplot(2,1,1)
[~,i] = max(x);
hold on
plot(x( (i-9):(i+10) ),'k');
plot(x( (i-9):(i+10) ),'k.');
box off
subplot(2,1,2)
[~,i] = max(y);
hold on
plot(y( (i-9):(i+10) ),'k');
plot(y( (i-9):(i+10) ),'k.');
box off

set(get(gcf,'Children'),'TickDir','out')
set(get(gcf,'Children'),'TickLength',[0.02 0.025])
set(get(gcf,'Children'),'YLim',[-4 4])
set(gcf,'Pos',[680   438   189   537])

%% autocorrelated
Npt = 1000;
Fs = 100;
x = randn(Npt,1);
y = randn(Npt,1);

mirrorpad = 0.5;
Whz = 20;
N = 3;

Wn = Whz / (Fs./2);
[b,a] = butter(N,Wn);
Npad = round(mirrorpad*Fs);

mirrordat = padarray(x, Npad, 'symmetric');
fltdat = filtfilt(b,a,mirrordat);
fltx = fltdat( (Npad+1):(end-Npad) );

mirrordat = padarray(y, Npad, 'symmetric');
fltdat = filtfilt(b,a,mirrordat);
flty = fltdat( (Npad+1):(end-Npad) );

%%
figure
subplot(2,1,1)
[~,i] = max(fltx);
hold on
plot(fltx( (i-9):(i+10) ),'k');
plot(fltx( (i-9):(i+10) ),'k.');
box off
subplot(2,1,2)
[~,i] = max(flty);
hold on
plot(flty( (i-9):(i+10) ),'k');
plot(flty( (i-9):(i+10) ),'k.');
box off

set(get(gcf,'Children'),'TickDir','out')
set(get(gcf,'Children'),'YLim',[-4 4])
set(gcf,'Pos',[680   438   189   537])
