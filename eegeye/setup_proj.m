

subjects = {'dgc','emb','gar','gvv','hab','hki','kah','r_k','s_h','swm','csj','drm',...
    'jfp','mip','stb'};

Ns = length(subjects);

fqnames={'1nc';'2c'}; % filtering conditions; 1nc - 1Hz non-causal, 2c - 2Hz causal

addpath('/home/robini/code/gauss_info')

data_dir = '/data/ssd/bubdetect';

time = -300:2:1000;
