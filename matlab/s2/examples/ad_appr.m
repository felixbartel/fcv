%
% this script calculates and minimizes the cv score for an quadrature
% example on the two-dimensional sphere
%
%% initialization

if ~exist('S2Fun','class') % add mtex library
  addpath '~/repo/mtex'; startup_mtex;
end
rng('default');                           % reset random generatormber of nodes
fun       = @(nodes) S2Fun.smiley(nodes); % example function
N         = 30;                           % bandwidth
nodes     = vector3d.rand(2*(N+1)^2);     % nodes in space domain
nodes     = unique(nodes(:));             % reshape to right format for nfsft
M         = length(nodes);                % number of nodes
W         = calcVoronoiArea(nodes);       % weights in space domain
f         = fun(nodes);                   % function values

f_e       = f+0.05*randn(size(f));        % noisy function values
s         = 3;                            % weights in requency domain
W_hat     = (2*(0:N)+1).^(2*s);
W_hat     = repelem(W_hat,1:2:(2*N+1))';


%% plot noisy data

pcolor(nodes,f_e,'contours',20,'upper'); hold on;
scatter(nodes,'MarkerSize',2,'MarkerColor','k'); hold off;
setColorRange([-0.5 0.5]);
mtexTitle('noisy data');
drawnow();


%% main computations

fcv = FCV_appr([nodes.rho,nodes.theta],f_e,W,N,s);

lambda_min = fcv.minimize();

[~,~,fhat_r] = fcv.compute(lambda_min);
sF_r = S2FunHarmonic(fhat_r);


%% plot reconstruction

nextAxis;
pcolor(sF_r,'contours',20,'upper');
setColorRange([-0.5 0.5]);
mtexTitle('reconstruction');
mtexColorbar;
