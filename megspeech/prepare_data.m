sub = 2;
data_dir = '/data1/ccni1/joachim/Projects/Story/Raw_cont';
ft_defaults

load(fullfile(data_dir,sprintf('Data_cont_%d_5.mat',sub)));
load('speech.mat')

%%
% Prepare MEG

% planar gradient
cfg = [];
cfg.method = 'distance';
cfg.neighbourdist = 0.037;
cfg.layout = '4D248';
cfg.feedback = 'no';
neighbours = ft_prepare_neighbours(cfg, data);

% FILTER
cfg = [];
cfg.padding = 0.2;
cfg.padtype = 'mirror';

cfg.lpfilter = 'yes';
cfg.lpfreq = 12;
cfg.lpfiltord = 3;
cfg.lpfilttype = 'but';

cfg.hpfilter = 'yes';
cfg.hpfreq = 2;
cfg.hpfiltord = 3;
cfg.hpfilttype = 'but';
    
flt_data = ft_preprocessing(cfg, data);

%%
    
% DOWNSAMPLE (full window to avoid edge artifacts)
cfg = [];
cfg.resamplefs = 50; % 20ms bins
cfg.detrend = 'no';
res_data = ft_resampledata(cfg, flt_data);

%%
% MEG PLANAR
cfg = [];
cfg.planarmethod = 'sincos';
cfg.neighbours = neighbours;
planar_data = ft_megplanar(cfg, res_data);

%%
% SAVE AS REGULAR ARRAY
plndat = cell2mat(reshape(planar_data.trial, [1 1 length(planar_data.trial)]));

%% filter speech
% low pass
mirrorpad = 0.2;
Whz = 12;
N = 3;
Fs = 250;

Wn = Whz / (Fs./2);
[b,a] = butter(N,Wn);

Npad = round(mirrorpad*Fs);
mirrordat = padarray(speech', Npad, 'symmetric');
fltspc = filtfilt(b,a,mirrordat);

% high pass filter 
Whz = 2;
N = 3;
Wn = Whz / (Fs./2);
[b,a] = butter(N,Wn,'high');
fltspc = filtfilt(b,a,fltspc);
fltspc = fltspc( (Npad+1):(end-Npad) );

fltspc = resample(fltspc, 1, 5);

%%
savefast sub2_dat plndat fltspc

%% convert channels to eeglab
chanlocs = [];
Nch = length(flt_data.label);
for ci=1:Nch
    thslab = flt_data.label{ci};
    chidx = find(strcmp(thslab,flt_data.grad.label));
    chanlocs(ci).labels = thslab;
    chanlocs(ci).X = flt_data.grad.chanpos(chidx,1);
    chanlocs(ci).Y = flt_data.grad.chanpos(chidx,2);
    chanlocs(ci).Z = flt_data.grad.chanpos(chidx,3); 
end

chanlocs = convertlocs(chanlocs,'cart2all');