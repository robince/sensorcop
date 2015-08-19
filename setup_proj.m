name = getComputerName();
if strcmp(name(1:10), 'robini2-pc')
    % robins workstation
    h.data_dir = '/data/hd/sensorcop';
else
    % on the grid
    h.data_dir = '/analyse/Project0110/robin/sensorcop';
end

addpath('/home/robini/code/gauss_info/')
addpath('/home/robini/code/gauss_info/mex/')
addpath('/home/robini/code/info/')
