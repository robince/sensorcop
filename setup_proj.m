name = getComputerName();
if strcmp(name(1:10), 'robini2-pc')
    % robins workstation
    data_dir = '/data/hd/sensorcop';
    addpath('/home/robini/code/gauss_info/')
    addpath('/home/robini/code/gauss_info/mex/')
    addpath('/home/robini/code/info/')
elseif strcmp(mexext,'mexmaci64')
    data_dir = '~/Documents/gladata/sensorcop';
    addpath('~/Documents/glacode/para_info/')
    addpath('~/Documents/glacode/para_info/mex/')
    addpath('~/Documents/glacode/info/')
else
    % on the grid
    data_dir = '/analyse/Project0110/robin/sensorcop';
    addpath('/home/robini/code/gauss_info/')
    addpath('/home/robini/code/gauss_info/mex/')
    addpath('/home/robini/code/info/')
end


