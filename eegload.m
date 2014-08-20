subj = ['HKI';'IJZ';'MAK';'S_H'];
Ns = length(subj);
filepath = '/analyse/Project0132/preprocessing/';
exptname='ft10';
fqnames = {'1nc';'2c'};

% for S = 1:Ns % Ns=4
S = 1;
        
%     for ses = 1:10
ses = 1;

subid = [subj(S,:) num2str(ses)];

FC = 1;

eeglab
filename=[exptname,'_',subid,'_ica_',fqnames{FC},'_clean.set'];
EEG = pop_loadset('filename',filename,'filepath',filepath);

% gather events
ev=[EEG.epoch.eventtype];nev=ev;
%fev=11:10:121;tev=12:10:122;
fev=11:10:121;tev=12:10:122;
for ee=1:length(ev)
    if ~isempty(intersect(ev(ee),fev))
        nev(ee)=1;
    else
        nev(ee)=2;
    end
end

% get data
ferp = EEG.csd_data(:,:,nev==1); % face ERPs
nerp = EEG.csd_data(:,:,nev==2); % noise ERPs

%%
dat = [];
dat.eegdat = EEG.data;
dat.csddat = EEG.csd_data;
dat.time = EEG.times;
dat.stim = nev-1; % 0 - face, 1 - noise
dat.chanlocs = EEG.chanlocs;

fname = sprintf('%s_ridat.mat',subid);
savefaststruct(fullfile(data_dir,fname),dat);