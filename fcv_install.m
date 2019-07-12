function fcv_install(manifold,nfft_path,mtex_path)
% FCV_INSTALL install the Fast Cros-validation toolbox
%
% Syntax:
%   s = fcv_install('t') for the torus
%   s = fcv_install('i') for the unit interval
%   s = fcv_install('s2') for the two-dmensional sphere
%   s = fcv_install('so3') for the rotation group
%   s = fcv_install('s2',nfft_path,mtex_path) if your library paths vary

if nargin < 3
  mtex_path = fullfile('~','repo','mtex');
end
if nargin < 2
  nfft_path = fullfile('~','repo','nfft');
end
if strcmp(manifold,'t')
  addpath(fullfile(nfft_path,'matlab','nfft'));
elseif strcmp(manifold,'i')
  addpath(fullfile(nfft_path,'matlab','nfct'));
elseif strcmp(manifold,'s2')
  addpath(fullfile(nfft_path,'matlab','nfsft'));
  addpath(mtex_path);
  startup_mtex
elseif strcmp(manifold,'so3')
  addpath(fullfile(nfft_dir,'matlab','nfsoft'));
  addpath(mtex_path);
  startup_mtex
end
folder = fileparts(which(mfilename('fullpath')));
addpath(genpath(fullfile(folder,manifold)))
end