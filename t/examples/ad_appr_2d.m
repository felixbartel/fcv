%
% this script calculates and minimizes the cv score for an nonequispaced
% example on the two-dimensional torus
%
%% initialization

rng('default');                        % reset random generator
fun       = @(x,y) peaks(6*x-3,6*y-3); % example function 
N         = 64;                        % bandwidth
nodes     = rand(2*N^2,2);             % nodes in space domain
nodes     = nodes.^2;
nodes     = uniquetol(nodes,1e-5,'ByRows',true);
f         = fun(nodes(:,1),nodes(:,2));% function values
f = f-min(f); f = f/max(f);

f_e       = f+0.05*randn(size(f));     % noisy function values
s         = 3;                         % weights in frequency domain


%% scatter noisy data

subplot(121);
scatter(nodes(:,1),-nodes(:,2),10,f_e,'filled');
title('noisy data');
axis square;
colorbar;
caxis([0 1]);
colormap jet;
drawnow();

%% main computations

fcv = FCV_appr(nodes,f_e,[],N,s);

lambda_min = fcv.minimize();

% calculate reconstruction
s = fcv.compute(lambda_min);

res = 480;
t = linspace(0,1,res);
[plotnodes_x,plotnodes_y] = meshgrid(t,t);
plan = nfft_init_2d(N,N,res^2);
nfft_set_x(plan,[plotnodes_x(:).';plotnodes_y(:).']);
nfft_precompute_psi(plan);
nfft_set_f_hat(plan,s.fhat_r);
nfft_trafo(plan);
plotf_r = nfft_get_f(plan);
nfft_finalize(plan);


%% plot reconstruction
subplot(122);
imagesc(real(reshape(plotf_r,res,res)));
title('reconstruction');
axis square;
colorbar;
caxis([0 1]);

