%
% this script calculates and plots the cv score for an quadrature
% example on the two-dimensional sphere
%
%% initialization

if ~exist('S2Fun','class') % add mtex library
  addpath '~/repo/mtex'; startup_mtex;
end
rng('default');                           % reset random generatormber of nodes
fun       = @(nodes) S2Fun.smiley(nodes); % example function
N         = 100;                          % bandwidth
[nodes,W] = quadratureS2Grid(2*N);        % nodes and weights in space domain
M         = length(nodes);                % number of nodes
f         = fun(nodes);                   % function values

f_e       = f+0.1*randn(size(f));         % noisy function values
s         = 3;                            % weights in requency domain
W_hat     = (2*(0:N)+1).^(2*s);
W_hat     = repelem(W_hat,1:2:(2*N+1))';

lambda    = 2.^(linspace(-38,-25,25));    % possible lambda
err       = 0*lambda;                     % stores L_2-error
cv        = 0*lambda;                     % stores cv score


%% main computations

nfsftmex('precompute',N,1000,1,0);
plan = nfsftmex('init_advanced',N,M,1);
nfsftmex('set_x',plan,[nodes.rho';nodes.theta']);
nfsftmex('precompute_x',plan);

f_hat = compute_f_hat(plan,f,0,W,W_hat);

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
    
% error
  f_hat_r = compute_f_hat(plan,f_e,lambda(idx),W,W_hat);
  err(idx) = norm(f_hat-f_hat_r);
  f_r = F(plan,f_hat_r);
  
% cv
  h_tilde = sum((2*(0:N)'+1)./(1+lambda(idx)*W_hat((1:N+1).^2)))/(4*pi)*W;
  cv(idx) = norm((f_r-f_e)./(1-h_tilde))^2;
end
close(wb);


%% computations for plotting

% get optimal lambda
[~,idx_err] = min(err);
[~,idx_cv]  = min(cv);

sF_r = S2FunHarmonic(compute_f_hat(plan,f_e,lambda(idx_cv),W,W_hat));

nfsftmex('finalize',plan);


%% plotting

figure(1);
% plot L_2-error
yyaxis left;
loglog(lambda,err); hold on;
scatter(lambda(idx_err),err(idx_err),40,'filled');
xlabel('\lambda');
ylabel('L_2-error');
% plot cv score
yyaxis right;
p = loglog(lambda,cv); hold on;
scatter(lambda(idx_cv),cv(idx_cv),40,'filled');
ylabel('cv score');
axis square;

% plot functions on the sphere
figure(2);
% plot noisy data
pcolor(nodes,f_e,'contours',20,'upper');
setColorRange([-0.5 0.5]);
%hold on; scatter(nodes,'MarkerSize',2,'MarkerColor','k'); hold off
mtexTitle('noisy data');
nextAxis;
% plot reconstruction
pcolor(sF_r,'contours',20,'upper');
setColorRange([-0.5 0.5]);
mtexTitle('reconstruction');
mtexColorbar;


%% helper functions

function f = F(plan,f_hat)
  nfsftmex('set_f_hat_linear',plan,f_hat);
  nfsftmex('trafo',plan);
  f = nfsftmex('get_f',plan);
  
  f = real(f);
end

function f_hat = compute_f_hat(plan,f,lambda,W,W_hat)
  f = W.*f;
  
  nfsftmex('set_f',plan, f);
  nfsftmex('adjoint',plan);
  f_hat = nfsftmex('get_f_hat_linear',plan);

  f_hat = f_hat./(1+lambda*W_hat);
end
