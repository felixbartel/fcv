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

lambda    = 2.^(linspace(-38,-25,25));    % possible lambda
err       = 0*lambda;                     % stores L_2-error
ocv       = 0*lambda;                     % stores ocv score
gcv       = 0*lambda;                     % stores gcv score
ocv_appr  = 0*lambda;                     % stores approximated ocv score
gcv_appr  = 0*lambda;                     % stores approximated gcv score

t_exact   = 0;
t_appr    = 0;


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
  tic;
  f_tilde = compute_f_hat(plan,f_e,lambda(idx),W,W_hat);
  nfsftmex('set_f_hat_linear',plan,f_tilde);
  nfsftmex('trafo',plan);
  f_r = nfsftmex('get_f',plan);
  f_r = real(f_r);
  t_tmp = toc;
  
  err(idx) = norm(f_hat-f_tilde);
  
% approximate cv
  tic;
  h_tilde = sum((2*(0:N)'+1)./(1+lambda(idx)*W_hat((1:N+1).^2)))/(4*pi)*W;
  ocv_appr(idx) = norm((f_r-f_e)./(1-h_tilde))^2;
  gcv_appr(idx) = norm((f_r-f_e)./(1-mean(h_tilde)))^2;
  t_appr = t_appr+t_tmp+toc;
  
% % exact cv
%   tic;
%   h = zeros(M,1);
%   for l = 1:M
%     tmp = double( 1:M == l )';
%     tmp = F(plan,compute_f_hat(plan,tmp,lambda(idx),W,W_hat));
%     h(l) = tmp(l);
%   end
%   ocv(idx) = norm((f_r-f_e)./(1-h))^2;
%   gcv(idx) = norm((f_r-f_e)./(1-mean(h)))^2;
%   t_exact = t_exact+t_tmp+toc;
end
close(wb);


%% computations for plotting

[~,idx_err]       = min(err);
[~,idx_gcv]       = min(gcv);
[~,idx_ocv]       = min(ocv);
[~,idx_gcv_appr]  = min(gcv_appr);
[~,idx_ocv_appr]  = min(ocv_appr);

sF_r = S2FunHarmonic(compute_f_hat(plan,f_e,lambda(idx_ocv_appr),W,W_hat));
nfsftmex('finalize',plan);

t_exact = t_exact/length(lambda);
t_appr  = t_appr/length(lambda);


%% plotting

figure(1);
% plot L_2-error
yyaxis left;
loglog(lambda,err); hold on;
scatter(lambda(idx_err),err(idx_err),'filled'); hold off;
xlabel('\lambda');
ylabel('L_2-error');
% plot cv score
yyaxis right;
p = loglog(lambda,[ocv; gcv; ocv_appr; gcv_appr]); hold on;
scatter(lambda([idx_ocv idx_gcv idx_ocv_appr idx_gcv_appr]),...
  [ocv(idx_ocv) gcv(idx_gcv) ocv_appr(idx_ocv_appr) gcv_appr(idx_gcv_appr)],40,'filled');
ylim([min([ocv gcv gcv_appr]) max([ocv gcv gcv_appr])]); hold off;
legend(p,'ocv','gcv','appr ocv','appr gcv');
ylabel('cv score');
axis square;

figure(2);
% plot noisy data
pcolor(nodes,f_e,'contours',20,'upper'); hold on;
scatter(nodes,'MarkerSize',2,'MarkerColor','k'); hold off;
setColorRange([-0.5 0.5]);
mtexTitle('noisy data');
nextAxis;
% plot reconstruction
pcolor(sF_r,'contours',20,'upper');
setColorRange([-0.5 0.5]);
mtexTitle('reconstruction');
mtexColorbar;


%% print times

fprintf('average time per lambda\n');
fprintf('\tt_exact = %f\n',t_exact);
fprintf('\tt_appr  = %f\n',t_appr);


%% helper functions

function f = F(plan,f_hat)
  nfsftmex('set_f_hat_linear',plan,f_hat);
  nfsftmex('trafo',plan);
  f = nfsftmex('get_f',plan);
  
  f = real(f);
end


function f_hat = compute_f_hat(plan,f,lambda,W,W_hat)
  f = sqrt(W).*f;
  f = [f;zeros(length(W_hat),1)];
  [f_hat,~] = lsqr(...
    @(x,transp_flag) afun(plan,x,lambda,W,W_hat,transp_flag),...
    f);
end

function y = afun(plan,x,lambda,W,W_hat,transp_flag)
  if strcmp(transp_flag,'notransp')
    nfsftmex('set_f_hat_linear',plan,x);
    nfsftmex('trafo',plan);
    y = nfsftmex('get_f',plan);
    
    y = sqrt(W).*y;
    
    y = [y; sqrt(lambda*W_hat).*x];
  else
    y = sqrt(W).*x(1:length(W));
    
    nfsftmex('set_f',plan,y);
    nfsftmex('adjoint',plan);
    y = nfsftmex('get_f_hat_linear',plan);
    
    y = y+sqrt(lambda*W_hat).*x(length(W)+1:end);
  end
end
