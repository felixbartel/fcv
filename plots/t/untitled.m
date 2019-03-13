%
% this script calculates and plots the cv score for an nonequispaced
% example on the one-dimensional torus
%
%% initialization

addpath '~/repo/nfft/matlab/nfft/'    % add the nfft-library
rng('default');                       % reset random generator
fun       = @(x) peaks(6*x-3,0);      % example function 
N         = 64;                       % bandwidth
nodes     = rand(2*N,1);              % nodes in space domain
nodes     = nodes.^2;
nodes     = unique(nodes,'rows');
M         = size(nodes,1);            % number of nodes
f         = fun(nodes);               % function values
f = f-min(f); f = f/max(f);           % normalize function

f_e       = f+0.05*randn(size(f));    % noisy function values
s         = 3;                        % weights in frequency domain
W_hat     = (abs((-N/2:N/2-1).').^s+1);

t         = linspace(-10,-3,50); % possible lambda
err       = 0*t;                 % stores L_2-error
ocv_appr  = 0*t;                 % stores approximated ocv score
gcv_appr  = 0*t;                 % stores approximated gcv score
grad_ocv_appr = 0*t;
grad_gcv_appr = 0*t;


%% original f_hat

M2 = 2^10;                        % sufficiently large bandwidth
[nodes2] = (0:M2-1)/M2;           % nodes in space domain
fe = fun(nodes2);                   % function values
fe = fe-min(fe); fe = fe/max(fe);
f_hat = ifft(fe);
f_hat = f_hat([(end-N/2+1):end 1:N/2]);


%% calculate Voronoi area (weights in space domain)

W = zeros(size(nodes));
W(1) = nodes(2)-nodes(end)+1;
W(2:end-1) = nodes(3:end)-nodes(1:end-2);
W(end) = nodes(1)+1-nodes(end-1);
W = W/2;


%% main computations

plan = nfft_init_1d(N,M);
nfft_set_x(plan,nodes');
nfft_precompute_psi(plan);

wb = waitbar(0);
for idx = 1:length(t) % loop over lambda
  waitbar(idx/length(t),wb);
  
  f_hat_r = compute_f_hat(plan,f_e,exp(t(idx)),W,W_hat);
  f_r = F(plan,f_hat_r);
  err(idx) = norm(f_hat-f_hat_r);
  
% approximated cv score
  h = sum(1./(1+exp(t(idx))*W_hat))*W;
  ocv_appr(idx) = log(norm((f_r-f_e)./(1-h))^2);
  gcv_appr(idx) = log(norm((f_r-f_e)./(1-mean(h)))^2);
  
  f_r_prime = F(-compute_f_hat(W_hat.*f_hat_r)); % lsqr multiplies with F\hermW^1/2 before what you don't want
  h_prime = -sum(W_hat./(1+exp(t(idx))*W_hat).^2)*W;
  
  grad_ocv_appr = 2*(((f_r-f).*f_r_prime)/(1-h).^2+(h_prime*(f_r-f).^2)/(1-h).^3);
  
end
close(wb);


%% computations for plotting

[~,idx_err]       = min(err);
[~,idx_gcv_appr]  = min(gcv_appr);
[~,idx_ocv_appr]  = min(ocv_appr);

f_hat_r = compute_f_hat(plan,f_e,exp(t(idx_ocv_appr)),W,W_hat);

nfft_finalize(plan);

res = 480;
plotnodes = linspace(0,1,res);

plan = nfft_init_1d(N,res);
nfft_set_x(plan,plotnodes);
nfft_precompute_psi(plan);
ploty_r = F(plan,f_hat_r);
nfft_finalize(plan);


%% plotting

subplot(121);
% plot noisy data
scatter(nodes,real(f_e),5,'k','filled'); hold on;
% plot reconstruction
plot(plotnodes,real(ploty_r)); hold off;
axis square;
title('noisy data and reconstruction');

subplot(222);
% plot L_2-error
yyaxis left;
loglog(exp(t),err); hold on;
scatter(exp(t(idx_err)),err(idx_err),'filled'); hold off;
xlabel('\lambda');
ylabel('L_2-error');
% plot cv score
yyaxis right;
p = loglog(exp(t),[ocv_appr; gcv_appr]); hold on;
scatter(exp(t([idx_ocv_appr idx_gcv_appr])),...
  [ocv_appr(idx_ocv_appr) gcv_appr(idx_gcv_appr)],40,'filled');
ylim([min([ocv_appr gcv_appr]) max(gcv_appr)]); hold off;
legend(p,'appr ocv','appr gcv');
ylabel('cv score');
axis square;

subplot(224);
axis square;

%% helper functions

function f = F(plan,f_hat)
  nfft_set_f_hat(plan,f_hat);
  nfft_trafo(plan);
  f = nfft_get_f(plan);
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
    nfft_set_f_hat(plan,x);
    nfft_trafo(plan);
    y = nfft_get_f(plan);
    
    y = sqrt(W).*y;
    
    y = [y; sqrt(lambda*W_hat).*x];
  else
    y = sqrt(W).*x(1:length(W));
    
    nfft_set_f(plan,y);
    nfft_adjoint(plan);
    y = nfft_get_f_hat(plan);
    
    y = y+sqrt(lambda*W_hat).*x(length(W)+1:end);
  end
end






