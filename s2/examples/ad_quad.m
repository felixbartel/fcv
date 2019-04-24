%
% this script calculates and minimizes the cv score for an quadrature
% example on the two-dimensional sphere
%
%% initialization

rng('default');                           % reset random generatormber of nodes
fun       = @(nodes) S2Fun.smiley(nodes); % example function
N         = 100;                          % bandwidth
[nodes,W] = quadratureS2Grid(2*N);        % nodes and weights in space domain
M         = length(nodes);                % number of nodes
f         = fun(nodes);                   % function values

f_e       = f+0.1*randn(size(f));         % noisy function values
s         = 3;                            % weights in requency domain


%% plot noisy data

pcolor(nodes,f_e,'contours',20,'upper');
setColorRange([-0.5 0.5]);
mtexTitle('noisy data');
drawnow();

%% main computations

fcv = FCV_quad([nodes.rho,nodes.theta],f_e,W,N,s);

lambda_min = fcv.minimize();
[~,~,fhat_r] = fcv.compute(lambda_min);
sF_r = S2FunHarmonic(fhat_r);


%% plot reconstruction

nextAxis;
pcolor(sF_r,'contours',20,'upper');
setColorRange([-0.5 0.5]);
mtexTitle('reconstruction');
mtexColorbar;
