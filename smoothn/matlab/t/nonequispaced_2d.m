%
% this script calculates the minimum of the cv score for an nonequispaced
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

lambda_0  = 1;
MaxFunEvals = 50;


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

[vertices,C] = voronoin(nodes_p);
W = zeros(M,1);
for idx = 1:M 
  W(idx) = polyarea(vertices(C{idx},1),vertices(C{idx},2));
end


%% scatter noisy data
clf; subplot(221);
scatter(nodes(:,1),-nodes(:,2),10,f_e,'filled');
title('noisy data');
axis square;
colorbar;
caxis(c);

colormap jet


%% main computations

plan = nfft_init_2d(N,N,M);
nfft_set_x(plan,nodes');
nfft_precompute_psi(plan);

options = optimset('OutputFcn',@outfun,'Display','iter','MaxFunEvals',MaxFunEvals);
subplot(122);
lambda_min = fminunc(@(lambda) V(plan, W, W_hat, f_e, lambda), ...
  lambda_0, options);
hold off;


%% computations for plotting

f_hat_r = compute_f_hat(plan,f_e,lambda_min,W,W_hat);

nfft_finalize(plan);

res = 480;
t = linspace(0,1,res);
[plotnodes_x,plotnodes_y] = meshgrid(t,t);
plan = nfft_init_2d(N,N,res^2);
nfft_set_x(plan,[plotnodes_x(:).';plotnodes_y(:).']);
nfft_precompute_psi(plan);
plotf_r = F(plan,f_hat_r);
nfft_finalize(plan);


%% plotting reconstruction

subplot(223);
imagesc(real(reshape(plotf_r,res,res)));
title('reconstruction');
axis square;
colorbar;
c = caxis;


%% helper functions

function stop = outfun(x, optimValues, ~)
  stop = false;
  loglog(x,optimValues.fval,'k.');
  xlabel('\lambda');
  ylabel('cv score');
  axis square;
  hold on;
  drawnow;
end

function ocv_appr = P(plan, W, W_hat, f, lambda)
  f_hat_r = compute_f_hat(plan,f,lambda,W,W_hat);
  f_r = F(plan,f_hat_r);
  
% approximated cv score
  h_tilde = sum(1./(1+lambda*W_hat))*W;
  ocv_appr = norm((f_r-f)./(1-h_tilde))^2;
end

function gcv_appr = V(plan, W, W_hat, f, lambda)
  f_hat_r = compute_f_hat(plan,f,lambda,W,W_hat);
  f_r = F(plan,f_hat_r);
  
% approximated cv score
  h_tilde = sum(1./(1+lambda*W_hat))*W;
  gcv_appr = norm((f_r-f)./(1-mean(h_tilde)))^2;
end

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
