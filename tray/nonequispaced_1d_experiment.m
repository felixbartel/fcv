%
% this script calculates and plots the cv score for an nonequispaced
% example on the one-dimensional torus
%
%% initialization

addpath '~/repo/nfft/matlab/nfft/'    % add the nfft-library
rng('default');                       % reset random generator
fun       = @(x) peaks(6*x-3,0);      % example function 
N         = 64;                       % bandwidth
nodes     = (1:(2*N))'./(2*N);

nodes = nodes(10:end);

nodes     = rand(2*N,1);              % nodes in space domain
nodes     = nodes.^2;
nodes     = unique(nodes,'rows');
M         = size(nodes,1);            % number of nodes
f         = fun(nodes);               % function values
f = f-min(f); f = f/max(f);           % normalize function

f_e       = f+0.05*randn(size(f));    % noisy function values
s         = 3;                        % weights in frequency domain
W_hat     = (abs((-N/2:N/2-1).').^s+1);

lambda    = 2.^(linspace(-16,-2,100)); % possible lambda
err       = 0*lambda;                 % stores L_2-error
ocv       = 0*lambda;                 % stores ocv score
gcv       = 0*lambda;                 % stores gcv score
ocv_appr  = 0*lambda;                 % stores approximated ocv score
gcv_appr  = 0*lambda;                 % stores approximated gcv score
ocv_cut   = 0*lambda;

hs = zeros(M,length(lambda));
hs_appr = hs;
hs_cut = hs;

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

%% exact quadrature
% 
% F_matrix = exp(2i*pi*nodes*(0:2*N));
% W = lsqr(F_matrix',double(0:2*N == 0)');
% W = 2*W;

%% main computations

plan = nfft_init_1d(N,M);
nfft_set_x(plan,nodes');
nfft_precompute_psi(plan);

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
  
  f_hat_r = compute_f_hat(plan,f_e,lambda(idx),W,W_hat);
  y_r = F(plan,f_hat_r);
  err(idx) = norm(f_hat-f_hat_r);
  
% approximated cv score
  h_tilde = sum(1./(1+lambda(idx)*W_hat))*W;
  ocv_appr(idx) = norm((y_r-f_e)./(1-h_tilde))^2;  
  gcv_appr(idx) = norm((y_r-f_e)./(1-mean(h_tilde)))^2;
  
  h_cut = h_tilde;
  h_cut( 0.9 < h_cut & h_cut <= 1 ) = 0.9;
  h_cut( 1 < h_cut & h_cut < 1.1 ) = 1.1;
  
  ocv_cut(idx) = norm((y_r-f_e)./(1-h_cut))^2;  
  
  hs_appr(:,idx) = h_tilde;  
  hs_cut(:,idx) = h_cut;
  
% exact cv score
  h = zeros(M,1);
  for l = 1:M
    tmp = double( 1:M == l )';
    tmp = F(plan,compute_f_hat(plan,tmp,lambda(idx),W,W_hat));
    h(l) = tmp(l);
  end
  ocv(idx) = norm((y_r-f_e)./(1-h))^2;
  gcv(idx) = norm((y_r-f_e)./(1-mean(h)))^2;
  
  hs(:,idx) = h;
  
end
close(wb);


%% computations for plotting

[~,idx_err]       = min(err);
[~,idx_gcv]       = min(gcv);
[~,idx_ocv]       = min(ocv);
[~,idx_gcv_appr]  = min(gcv_appr);
[~,idx_ocv_appr]  = min(ocv_appr);

f_hat_r = compute_f_hat(plan,f_e,lambda(idx_ocv_appr),W,W_hat);

nfft_finalize(plan);

res = 480;
plotnodes = linspace(0,1,res);

plan = nfft_init_1d(N,res);
nfft_set_x(plan,plotnodes);
nfft_precompute_psi(plan);
ploty_r = F(plan,f_hat_r);
nfft_finalize(plan);


%% plotting

clf;
subplot(221);
% plot noisy data
scatter(nodes,real(f_e),5,'k','filled'); hold on;
% plot reconstruction
plot(plotnodes,real(ploty_r)); hold off;
%axis square;
title('noisy data and reconstruction');

subplot(222);
% plot L_2-error
yyaxis left;
loglog(lambda,err); hold on;
%scatter(lambda(idx_err),err(idx_err),'filled'); hold off;
xlabel('\lambda');
ylabel('L_2-error');
% plot cv score
yyaxis right;
p = loglog(lambda,[ocv; ocv_appr; ocv_cut]); hold on;
% p = loglog(lambda,[ocv; gcv; ocv_appr; gcv_appr]); hold on;
% scatter(lambda([idx_ocv idx_gcv idx_ocv_appr idx_gcv_appr]),...
%   [ocv(idx_ocv) gcv(idx_gcv) ocv_appr(idx_ocv_appr) gcv_appr(idx_gcv_appr)],40,'filled');
% ylim([min([ocv gcv gcv_appr]) max([ocv gcv gcv_appr])]); hold off;
% legend(p,'ocv','gcv','appr ocv','appr gcv');
legend(p,'ocv','ocv appr','ocv cut');
ylabel('cv score');
%axis square;

subplot(224);
semilogx(lambda,0*lambda+1,'k'); hold on;
tmp = ones(M,1)*lambda;
scatter(tmp(:),hs(:),10,'b','filled','MarkerFaceAlpha',0.1);
scatter(tmp(:),hs_cut(:),10,'r','filled','MarkerFaceAlpha',0.1);
scatter(tmp(:),hs_appr(:),10,'k','filled','MarkerFaceAlpha',0.1); hold off;
%axis square;

subplot(223);
F_matrix = exp(2i*pi*nodes*(0:N));
pcolor(abs(F_matrix'*diag(W)*F_matrix-eye(N+1)));
colorbar;


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






