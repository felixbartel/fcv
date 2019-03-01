%
% this script calculates and plots the cv score for an nonequispaced
% example on the two-dimensional torus
%
%% initialization

addpath '~/repo/nfft/matlab/nfft/'     % add the nfft-library
rng('default');                        % reset random generator
fun       = @(x,y) peaks(6*x-3,6*y-3); % example function 
N         = 64;                        % bandwidth
nodes     = rand(2*N^2,2);             % nodes in space domain
nodes     = nodes.^2;
nodes     = uniquetol(nodes,1e-5,'ByRows',true);
M         = size(nodes,1);             % number of nodes 

f         = fun(nodes(:,1),nodes(:,2));% function values
f = f-min(f); f = f/max(f);

f_e       = f+0.05*randn(size(f));     % noisy function values
s         = 3;                         % weights in frequency domain
W_hat = (1+(-N/2:N/2-1)'.^2+(-N/2:N/2-1).^2).^(s/2);
W_hat = W_hat(:);

lambda    = 2.^(linspace(-12,-6,25));  % possible lambda
err       = 0*lambda;                  % stores L_2-error
ocv       = 0*lambda;                  % stores ocv score
gcv       = 0*lambda;                  % stores approximated ocv score
ocv_appr  = 0*lambda;                  % stores gcv score
gcv_appr  = 0*lambda;                  % stores approximated gcv score

t_exact   = 0;
t_appr    = 0;


%% original f_hat

M2 = 2^10;                                              % sufficiently large bandwidth
[nodes_x,nodes_y] = meshgrid((0:M2-1)/M2,(0:M2-1)/M2);  % nodes in space domain
fe = reshape(fun(nodes_x,nodes_y),[],1);                  % function values
fe = fe-min(fe); fe = fe/max(fe);

f_hat = 1/M2^2*fft2(fe,'transp');
f_hat = reshape(f_hat,M2,M2);
f_hat = f_hat([(M2-N/2+1):M2 1:(N/2)],[(M2-N/2+1):M2 1:(N/2)]);
f_hat = f_hat(:);


%% calculate Voronoi area
nodes_p = [nodes; ...
  nodes(:,1)-1 nodes(:,2)-1; ...
  nodes(:,1)-1 nodes(:,2)+0; ...
  nodes(:,1)-1 nodes(:,2)+1; ...
  nodes(:,1)+0 nodes(:,2)-1; ...
  nodes(:,1)+0 nodes(:,2)+1; ...
  nodes(:,1)+1 nodes(:,2)-1; ...
  nodes(:,1)+1 nodes(:,2)+0; ...
  nodes(:,1)+1 nodes(:,2)+1];

[V,C] = voronoin(nodes_p);
W = zeros(M,1);
for idx = 1:M 
  W(idx) = polyarea(V(C{idx},1),V(C{idx},2));
end

% % draw Voronoi decomposition
% for j = 1:N
%   patch(V(C{j},1),V(C{j},2),W(j));
% end
% axis([0 1 0 1]);
% axis square;


%% main computations

plan = nfft_init_2d(N,N,M);
nfft_set_x(plan,nodes');
nfft_precompute_psi(plan);

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
  tic;
  f_hat_r = compute_f_hat(plan,f_e,lambda(idx),W,W_hat);
  f_r = F(plan,f_hat_r);
  t_tmp = toc;
  err(idx) = norm(f_hat-f_hat_r);
  
% approximated cv score
  tic;
  h_tilde = sum(1./(1+lambda(idx)*W_hat))*W;
  ocv_appr(idx) = norm((f_r-f_e)./(1-h_tilde))^2;  
  gcv_appr(idx) = norm((f_r-f_e)./(1-mean(h_tilde)))^2;
  t_appr = t_appr+t_tmp+toc;
  
% % exact cv score
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

f_hat_r = compute_f_hat(plan,f_e,lambda(idx_ocv_appr),W,W_hat);

nfft_finalize(plan);

res = 480;
t = linspace(0,1,res);
[plotnodes_x,plotnodes_y] = meshgrid(t,t);
plan = nfft_init_2d(N,N,res^2);
nfft_set_x(plan,[plotnodes_x(:).';plotnodes_y(:).']);
nfft_precompute_psi(plan);
plotf_r = F(plan,f_hat_r);
nfft_finalize(plan);


%% plotting

% plot reconstruction
subplot(223);
imagesc(real(reshape(plotf_r,res,res)));
title('reconstruction');
axis square;
colorbar;
c = caxis;

% scatter noisy data
subplot(221);
scatter(nodes(:,1),-nodes(:,2),10,f_e,'filled');
title('noisy data');
axis square;
colorbar;
caxis(c);

colormap jet

subplot(122);
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


%% print times

fprintf('average time per lambda\n');
fprintf('\tt_exact = %f\n',t_exact/length(lambda));
fprintf('\tt_appr  = %f\n',t_appr/length(lambda));


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
