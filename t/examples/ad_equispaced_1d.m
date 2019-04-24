%
% this script minimizes the cv score for an equispaced example on the one-
% dimensional torus
%
%% initialization

rng('default');                      % reset random generator

fun    = @(x) peaks(6*x-3,0);        % example function
M      = 2^10;                       % number of nodes
nodes  = (0:M-1)'/M;                 % nodes in space domain
f      = fun(nodes);                 % function values
f = f-min(f); f = f/max(f);          % normalize function

f_e    = f+0.1*randn(size(f));       % noisy function values
s      = 3;                          % weights in frequency domain


%% plot noisy data

scatter(nodes,real(f_e),5,'k','filled'); hold on;
axis square;
title('noisy data');
drawnow();


%% main computations

fcv = FCV_equispaced(1,f_e,s);
lambda_min = fcv.minimize();

% calculate reconstruction
[~,~,f_r] = fcv.compute(lambda_min);


%% plot reconstruction

plot(nodes,real(f_r)); hold off;
title('noisy data and reconstruction');