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
nu        = 3;
W_hat     = ((1:M).').^nu;                  % weights in frequency domain

lambda    = 2.^linspace(-17,-11,25);        % possible lambda
err       = 0*lambda;                       % stores l_2-error
ocv       = 0*lambda;                       % stores ocv score
gcv       = 0*lambda;                       % stores gcv score
ocv_appr  = 0*lambda;                       % stores approximated ocv score
gcv_appr  = 0*lambda;                       % stores approximated gcv score

t_exact   = 0;
t_appr    = 0;

%% calculate weights in space domain

[nodes_tilde,idx] = sort(nodes);
nodes_tilde = acos(nodes_tilde);

W = pi-(nodes_tilde(1)+nodes_tilde(2))/2;               % first node
W = [W; (nodes_tilde(1:end-2)-nodes_tilde(3:end))/2];   % indermediate nodes
W = [W; (nodes_tilde(end-1)+nodes_tilde(end))/2];       % last nodes
W(idx) = 2*W/pi;

% scatter(nodes, zeros(size(nodes)), 50, W, 'filled')
% colorbar;


%% main computations

plan = nfct_init_1d(M,length(nodes));
nfct_set_x(plan,acos(nodes.')/(2*pi));

plan2 = nfct_init_1d(2*M,length(nodes));
nfct_set_x(plan2,acos(nodes.')/(2*pi));

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
  
  tic;
  f_r = mult_H(plan,f_e,lambda(idx),W,W_hat);
  err(idx) = norm(f-f_r);
  t_tmp = toc;
  
% compute diagonal emelents
  tic;
  a = 1./(1+lambda(idx)*W_hat);
  a(1) = a(1)/2;

  a_hat_r = zeros(2*M,1);
  a_hat_r(1:2:2*M-1) = a;

  nfct_set_f_hat(plan2,double(a_hat_r));
  nfct_trafo(plan2);
  h = nfct_get_f(plan2);

  h = W.*(h+sum(a))/2;

  ocv_appr(idx) = norm((f_r-f_e)./(1-h))^2;
  gcv_appr(idx) = norm((f_r-f_e)./(1-mean(h)))^2;
  t_appr = t_appr+t_tmp+toc;
  
% % exact cv score
%   tic;
%   h = zeros(M,1);
%   for l = 1:M
%     tmp = double( 1:M == l )';
%     tmp = mult_H(plan,tmp,lambda(idx),W,W_hat);
%     h(l) = tmp(l);
%   end
%   ocv(idx) = norm((f_r-f_e)./(1-h))^2;
%   gcv(idx) = norm((f_r-f_e)./(1-mean(h)))^2;
%   t_exact = t_exact+t_tmp+toc;
end
close(wb);


%% computations for plotting

[~,idx_gcv] = min(gcv);
[~,idx_ocv] = min(ocv);
[~,idx_gcv_appr] = min(gcv_appr);
[~,idx_ocv_appr] = min(ocv_appr);

f_hat_r = sqrt(M/2)*sqrt(W).*f_e;
f_hat_r = [f_hat_r;zeros(length(W_hat),1)];
[f_hat_r,~] = lsqr(...
  @(x,transp_flag) afun(plan,x,lambda(idx_ocv_appr),W,W_hat,transp_flag),...
  f_hat_r,1e-8,1000);

nfct_finalize(plan);
nfct_finalize(plan2);

t_exact = t_exact/length(lambda);
t_appr = t_appr/length(lambda);

% calculate values for plotting
res = 480;
plotnodes = linspace(-1,1,res);

plan = nfct_init_1d(M,length(plotnodes));
nfct_set_x(plan,acos(plotnodes)/(2*pi));

plotf_r = ndctIII(plan,f_hat_r);  
nfct_finalize(plan);


%% plotting

[nodes,idx] = sort(nodes);
subplot(121);
% plot noisy data
scatter(nodes,real(f_e(idx)),10,'filled','k'); hold on;
% plot reconstruction
plot(plotnodes,real(plotf_r)); hold off;
axis square;
hold off;
title('noisy data and reconstruction');

subplot(122);
% plot l_2-error
yyaxis left;
loglog(lambda,err);
xlabel('\lambda');
ylabel('l_2-error');
% plot cv score
yyaxis right;
p = loglog(lambda,[ocv;gcv;ocv_appr;gcv_appr]); hold on;
scatter(lambda([idx_ocv idx_gcv idx_ocv_appr idx_gcv_appr]),...
  [ocv(idx_ocv) gcv(idx_gcv) ocv_appr(idx_ocv_appr) gcv_appr(idx_gcv_appr)],40,'filled'); hold off;
ylim([min([ocv gcv gcv_appr]) max([ocv gcv gcv_appr])]); hold off;
legend(p,'ocv','gcv','appr ocv','appr gcv');
ylabel('cv score');
axis square;
xlim([lambda(1) lambda(end)]);


%% print times

fprintf('\naverage time per lambda\n');
fprintf('\tt_exact = %f\n',t_exact);
fprintf('\tt_appr  = %f\n',t_appr);


%% helper functions

function f_hat = mult_H(plan,f,lambda,W,W_hat)
  M = length(f);
  f = sqrt(M/2*W).*f;
  f = [f;zeros(length(W_hat),1)];
  [f_hat,~] = lsqr(...
    @(x,transp_flag) afun(plan,x,lambda,W,W_hat,transp_flag),...
    f,1e-8,1000);
  
  f_hat = ndctIII(plan,f_hat); 
end

function y = afun(plan,x,lambda,W,W_hat,transp_flag)
  M = length(W);
  if strcmp(transp_flag,'notransp')
    y = sqrt(M/2*W).*ndctIII(plan,x);
    y = [y; sqrt(lambda*W_hat).*x];
  else
    y = ndctII(plan,sqrt(M/2*W).*x(1:M));
    y = y+sqrt(lambda*W_hat).*x(M+1:end);
  end
end

function f_hat = ndctII(plan,f)
  M = length(f);
  nfct_set_f(plan,double(f));
  nfct_adjoint(plan);
  f_hat_adjoint = nfct_get_f_hat(plan);

  f_hat = sqrt(2/M)*f_hat_adjoint;
  f_hat(1) = f_hat(1)/sqrt(2);
end

function f = ndctIII(plan,f_hat)
  M = length(f_hat);
  f_hat = [f_hat(1)/sqrt(2); f_hat(2:end)];
  nfct_set_f_hat(plan,double(f_hat));
  nfct_trafo(plan);
  f = nfct_get_f(plan);

  f = sqrt(2/M)*f;
end
