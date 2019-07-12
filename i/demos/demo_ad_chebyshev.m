%
% this script calculates and plots the cv score for an example on the unit
% interval with chebyshev nodes
%
%% initialization

rng('default');                       % reset random generator
fun    = @(x) peaks(2*x,0);           % example function
M      = 2^7;                         % number of nodes
nodes  = cos(pi*(2*(1:M)-1)/(2*M)).'; % Tschebyshev nodes of first kind
f      = fun(nodes);                  % function values
f = f-min(f); f = f/max(f);           % normalize function

f_e    = f+0.05*randn(size(f));       % noisy function values
s      = 3;


%% plot noisy data

scatter(nodes,real(f_e),10,'k','filled'); hold on;
axis square;
title('noisy data');
drawnow();


%% main computations

fcv = FCV_chebyshev(f_e,s);
lambda_min = fcv.minimize();
res = fcv.compute(lambda_min);


%% plot reconstruction

plot(nodes,real(res.f_r)); hold off;
title('noisy data and recontruction')
