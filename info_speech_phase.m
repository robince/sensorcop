
% load data

dat = load(fullfile(data_dir,'forRobin.mat'));
Fs = dat.sf;

spc = dat.dat(1,:);
meg = dat.dat(2,:);

%% band pass filter both
mirrorpad = 1;
Whz = [4 10];
N = 3;

Wn = Whz / (Fs./2);
[b,a] = butter(N,Wn);

Npad = round(mirrorpad*Fs);

mirrordat = padarray(spc', Npad, 'symmetric');
fltdat = filtfilt(b,a,mirrordat);
fltspc = fltdat( (Npad+1):(end-Npad) );
hspc = hilbert(fltspc);

mirrordat = padarray(meg', Npad, 'symmetric');
fltdat = filtfilt(b,a,mirrordat);
fltmeg = fltdat( (Npad+1):(end-Npad) );
hmeg = hilbert(fltmeg);

optlag = 34; % from time domain cross-info
hspc = hspc(1:(end-optlag));
fltspc = fltspc(1:(end-optlag));
dfltspc = gradient(fltspc);

hmeg = hmeg((optlag+1):end);
fltmeg = fltmeg((optlag+1):end);
dfltmeg = gradient(fltmeg); 

%%
cfltspc = copnorm(fltspc);
cfltmeg = copnorm(fltmeg);
cdfltspc = copnorm(dfltspc);
cdfltmeg = copnorm(dfltmeg);

% Iraw = info_gg(cfltspc, cfltmeg, true, false, false);

zcopnorm = @(x) copnorm([real(x) imag(x)]);
chmeg = zcopnorm(hmeg);
chspc = zcopnorm(hspc);

phsspc = (hspc ./ abs(hspc));% .* rand(size(hspc));
cphsspc = zcopnorm(phsspc);
% clinphsspc = copnorm(angle(hspc));
% powspc = abs(hspc) .* exp(i.*( (rand(size(hspc))*2*pi) - pi ));
% cpowspc = zcopnorm(powspc);
cpowspc = copnorm(abs(hspc));

phsmeg = (hmeg ./ abs(hmeg));% .* rand(size(hmeg));
cphsmeg = zcopnorm(phsmeg);
% clinphsmeg = copnorm(angle(hmeg));
% powmeg = abs(hmeg) .* exp(i.*( (rand(size(hmeg))*2*pi) - pi ));
% cpowmeg = zcopnorm(powmeg);
cpowmeg = copnorm(abs(hmeg));

Iphs = info_gg(cphsspc, cphsmeg, true, false, false)
% Ipow = info_gg(cpowspc, cpowmeg, true, false, false)
% Ih = info_gg(chspc, chmeg, true, false, false)
% Ipp = info_gg([cpowspc cphsspc], [cpowmeg cphsmeg], true, false, false)

Ihphs = info_gg(chspc, cphsmeg, true, false, false)
Ippphs = info_gg([cpowspc cphsspc], cphsmeg, true, false, false)
Ipowphs = info_gg([cpowspc], cphsmeg, true, false, false)

Iflt_phs = info_gg(cfltspc, cphsmeg, true,false,false)
Ifltphs_phs = info_gg([cfltspc cphsspc], cphsmeg, true,false,false)
Ifltpp_phs = info_gg([cfltspc cpowspc cphsspc], cphsmeg, true,false,false)
Iflth_phs = info_gg([cfltspc chspc], cphsmeg, true,false,false)
% Ilinphs = info_gg(clinphsspc, clinphsmeg, true, false, false)
% Ihlinphase = info_gg(chspc, clinphsmeg, true, false, false)
% Iplplinphase = info_gg([clinphsspc cpowspc], clinphsmeg, true, false, false)
% 
% Ilinphsh = info_gg(clinphsspc, clinphsmeg, true, false, false)
% Ihlinphase = info_gg(chspc, clinphsmeg, true, false, false)
% Iplplinphase = info_gg([clinphsspc cpowspc], clinphsmeg, true, false, false)

%%


Ipowphs = info_gg(cpowspc, cphsmeg, true, false, false)
Ipowh = info_gg(cpowspc, chmeg, true, false, false)

Iphspow = info_gg(cphsspc, cpowmeg, true, false, false)
Iphsh = info_gg(cphsspc, chmeg, true, false, false)





% Irawphs = info_gg(cfltspc, cphsmeg, true, false, false)
% Irawpow = info_gg(cfltspc, cpowmeg, true, false, false)
% Irawh = info_gg(cfltspc, chmeg, true, false, false)

%% bootstrap phase
% Nboot = 1000;
% Iboot = zeros(1,Nboot);
% 
% parfor bi=1:Nboot
%     Iboot(bi) = info_gg(cphsspc(randperm(size(cphsspc,1)),:), cphsmeg, true, false, false);
% %     Iboot(bi) = info_gg(chspc(randperm(size(cphsspc,1)),:), chmeg, true, false, false);
% end
% %%
% max(Iboot)

%%
ccfltspc = normcdf(cfltspc);
ccfltmeg = normcdf(cfltmeg);
Nplt = 50000;
pltidx = randperm(length(cfltspc),Nplt);
figure
scatter(ccfltspc(pltidx),ccfltmeg(pltidx),6,'k', 'filled')
