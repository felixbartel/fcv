%
% this script minimizes the cv score for an equispaced example on the two-
% dimensional torus
%
%% initialization

rng('default');                                     % reset random generator

fun    = @(x,y) peaks(6*x-3,6*y-3);                 % example function
N      = 2^10;                                      % bandwidth
[nodes_x,nodes_y] = meshgrid((0:N-1)/N,(0:N-1)/N);  % nodes in space domain
f      = reshape(fun(nodes_x,nodes_y),[],1);        % function values
f = f-min(f); f = f/max(f);                         % normalize function
fhat   = fftd_adj(f,2)/N^2;                         % get original fhat

f_e    = f+0.1*randn(size(f));                      % noisy function values
s      = 3;                                         % weights in frequency domain


%% plot noisy data

subplot(121);
imagesc(real(reshape(f_e,N,N)));
title('noisy data');
axis square;
caxis([0 1]);
colorbar;
colormap jet;
drawnow();


%% main computations

fcv = FCV_equispaced(2,f_e,s);
lambda_min = fcv.minimize();

% calculate reconstruction
s = fcv.compute(lambda_min);


%% plot reconstruction

subplot(122);
imagesc(real(reshape(s.f_r,N,N)));
title('reconstruction');
axis square;
caxis([0 1]);
colorbar;