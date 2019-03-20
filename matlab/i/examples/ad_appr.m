%
% this script calculates and plots the cv score for an example on the unit
% interval with nonequispaced nodes
%
%% initialization


addpath('~/repo/nfft/matlab/nfct/');        % add ndct library
rng('default');                             % reset random generator
fun       = @(x) peaks(2*x,0);              % example function
M         = 2^7;                            % number of nodes
nodes     = linspace(-1,1,M)';              % nodes in space domain
nodes     = nodes+0.01*randn(size(nodes));
if nnz(( -1 > nodes ) | ( nodes > 1 ))
  nodes(( -1 > nodes ) | ( nodes > 1 )) = 2*rand(nnz(( -1 > nodes ) | ( nodes > 1 )),1)-1;
end
f         = fun(nodes);                     % function values
f = f-min(f); f = f/max(f);                 % normalize function

f_e       = f+0.05*randn(size(f));          % noisy function values
s         = 3;


%% plot noisy data

[~,idx] = sort(nodes);
scatter(nodes(idx),real(f_e(idx)),10,'filled','k'); hold on;
axis square;
title('noisy data');
drawnow();


%% main computations

fcv = FCV_appr(nodes,f_e,[],M,s);
lambda_min = fcv.minimize();
[~,~,fhat_r] = fcv.compute(lambda_min);

res = 480;
plotnodes = linspace(-1,1,res);

plan = nfct_init_1d(M,length(plotnodes));
nfct_set_x(plan,acos(plotnodes)/(2*pi));

plotf_r = sqrt(M/2)*ndctIII(plan,fhat_r);
nfct_finalize(plan);


%% plot reconstruction

plot(plotnodes,real(plotf_r)); hold off;
title('noisy data and reconstruction');
