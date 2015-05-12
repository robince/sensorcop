% load data

dat = load(fullfile(data_dir,'forRobin.mat'));
Fs = dat.sf;

spc = dat.dat(1,:);
meg = dat.dat(2,:);

%% low pass filter both
mirrorpad = 0.5;
Whz = 20;
N = 3;

Wn = Whz / (Fs./2);
[b,a] = butter(N,Wn);

Npad = round(mirrorpad*Fs);

mirrordat = padarray(spc', Npad, 'symmetric');
fltdat = filtfilt(b,a,mirrordat);
fltspc = fltdat( (Npad+1):(end-Npad) );

mirrordat = padarray(meg', Npad, 'symmetric');
fltdat = filtfilt(b,a,mirrordat);
fltmeglp = fltdat( (Npad+1):(end-Npad) );

%% high pass filter MEG
Whz = 0.8;
N = 3;
Wn = Whz / (Fs./2);
[b,a] = butter(N,Wn,'high');
mirrordat = padarray(fltmeglp, Npad, 'symmetric');
fltdat = filtfilt(b,a,mirrordat);
fltmeg = fltdat( (Npad+1):(end-Npad) );

%% simple x-correlation
[xc,lags] = xcorr(fltspc, fltmeg, 100);
% plot(lags./Fs, xc);

Nlags = length(lags);
% dfltspc = gradient(fltspc);
% dfltmeg = gradient(fltmeg);
% qspc = copnorm(dfltspc);
% qmeg = copnorm(dfltmeg);

qspc = copnorm(fltspc);
qmeg = copnorm(fltmeg);

Ic = zeros(Nlags,1);
for li=1:Nlags
    lag = lags(li);
    if lag==0
        x = qspc;
        y = qmeg;
    elseif lag<0 % speech leads
        l = abs(lag);
        x = qspc(1:(end-l));
        y = qmeg((l+1):end);
    elseif lag>0 % meg leads
        x = qspc((l+1):end);
        y = qmeg(1:(end-l));
    end
    Ic(li) = info_gg(x,y,true,true,false);
end

%%
figure
subplot(211)
plot(lags./Fs, xc);
title('XCORR')
subplot(212)
plot(lags./Fs, Ic)
title('Info')

%%
% max effect
lag = -34;



