addpath '~/repo/nfft/matlab/nfft'
folder = fileparts(which(mfilename('fullpath')));
addpath(folder(1:end-2))
addpath(fullfile(folder,'examples'))
addpath(fullfile(folder,'tools'))
