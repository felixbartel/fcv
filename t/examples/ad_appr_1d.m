%
% this script calculates and minimizes the cv score for an nonequispaced
% example on the one-dimensional torus
%
%% initialization

rng('default');                       % reset random generator
fun       = @(x) peaks(6*x-3,0);      % example function 
N         = 2^6;                     % bandwidth
nodes     = rand(2*N,1);              % nodes in space domain
nodes     = nodes.^2;
nodes     = unique(nodes,'rows');
f         = fun(nodes);               % function values
f = f-min(f); f = f/max(f);           % normalize function

f_e       = f+0.05*randn(size(f));    % noisy function values
s         = 3;                        % weights in frequency domain


%% plot noisy data

scatter(nodes,real(f_e),5,'k','filled'); hold on;
axis square;
title('noisy data');
drawnow();


%% main computations

fcv = FCV_appr(nodes,f_e,[],N,s);
lambda_min = fcv.minimize();

% calculate reconstruction
[~,~,fhat_r] = fcv.compute(lambda_min);

res = 480;
plotnodes = linspace(0,1,res);

plan = nfft_init_1d(N,res);
nfft_set_x(plan,plotnodes);
nfft_precompute_psi(plan);
nfft_set_f_hat(plan,fhat_r);
nfft_trafo(plan);
plotf_r = nfft_get_f(plan);
nfft_finalize(plan);

%% plot reconstruction

plot(plotnodes,real(plotf_r)); hold off;
axis square;
title('noisy data and reconstruction');