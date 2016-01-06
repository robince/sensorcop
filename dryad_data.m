dryad_dir = '/home/robini/elife_dryad/data';

dat = load(fullfile(dryad_dir,'eeg_face_vs_noise'));

time = dat.time;
idx = (time>-200) & (time<600);
dat.csddat = single(dat.csddat(:,idx,:));
dat.time = time(idx);

%%
savefaststruct(fullfile(dryad_dir,'eeg_face_vs_noise'), dat)


%%
dryad_dir = '/home/robini/elife_dryad/data';
datfull = load(fullfile(dryad_dir,'data_dgc_1nc'));

dat = [];
idx = (time>-200) & (time<600);
dat.time = time(idx);
dat.csddat = datfull.ferp(:,idx,datfull.REmaxMI);
dat.chanlocs = datfull.chanlocs;
dat.electrode = datfull.REmaxMI;
dat.eyestim = datfull.eyebubs;

%%
savefaststruct(fullfile(dryad_dir,'eeg_eye_visibility'), dat)