%
% this script calculates and plots the cv score for an nonequispaced
% example on the one-dimensional torus
%
%% initialization

addpath '~/repo/nfft/matlab/nfft/'    % add the nfft-library
addpath '~/repo/nfft/matlab/infft1d/' % add the infft1d-library
rng('default');                       % reset random generator
fun       = @(x) peaks(6*x-3,0);      % example function 
N        = 64;                       % bandwidth
nodes     = rand(2*N,1);             % nodes in space domain
nodes     = nodes.^2;
nodes     = unique(nodes,'rows');

M         = 2*N;                     % number of nodes  
% Jittered nodes in [-0.5,0.5)
nodes = (-0.5:1/M:0.5-1/M) + 1/(4*M)*rand(1,M);
nodes = nodes.'+0.5;
                
f         = fun(nodes);                 % function values
f = f-min(f); f = f/max(f);           % normalize function

f_e       = f+0.05*randn(size(f));    % noisy function values
s         = 3;                        % weights in frequency domain
W_hat     = (abs((-N/2:N/2-1).').^s+1);

lambda    = 2.^(linspace(-20,-5,25)); % possible lambda
err       = 0*lambda;                 % stores L_2-error
ocv       = 0*lambda;                 % stores ocv score
gcv       = 0*lambda;                 % stores gcv score
ocv_appr  = 0*lambda;                 % stores approximated ocv score
gcv_appr  = 0*lambda;                 % stores approximated gcv score

t_exact   = 0;
t_appr    = 0;


%% original f_hat

M2 = 2^10;                        % sufficiently large bandwidth
[nodes2] = (0:M2-1)/M2;           % nodes in space domain
fe = fun(nodes2);                   % function values
fe = fe-min(fe); fe = fe/max(fe);
f_hat = ifft(fe);
f_hat = f_hat([(end-N/2+1):end 1:N/2]).';


%% main computations

nodes_nfft = mod(nodes'-0.5,1)-0.5;

plan_nfft = nfft_init_1d(N,M);
nfft_set_x(plan_nfft,nodes_nfft);
nfft_precompute_psi(plan_nfft);

plan_infft = infft(-nodes_nfft,N);

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
  
  tic;
  f_hat_r = compute_f_hat(plan_nfft,plan_infft,f_e,lambda(idx),W_hat);
  f_r = F(plan_nfft,f_hat_r);
  t_tmp = toc;
  err(idx) = norm(f_hat-f_hat_r);
  
% approximated cv score
  tic;
  h = zeros(M,1);
  for l = 1:M
    tmp = double( 1:M == l )';
    tmp = F(plan_nfft,compute_f_hat_appr(plan_nfft,plan_infft,tmp,lambda(idx),W_hat));
    h(l) = tmp(l);
  end
  ocv_appr(idx) = norm((f_r-f_e)./(1-h))^2;
  gcv_appr(idx) = norm((f_r-f_e)./(1-mean(h)))^2;
  t_appr = t_appr+t_tmp+toc;
%   
% % exact cv score
%   tic;
%   h = zeros(M,1);
%   for l = 1:M
%     tmp = double( 1:M == l )';
%     tmp = F(plan_nfft,compute_f_hat(plan_nfft,plan_infft,tmp,lambda(idx),W_hat));
%     h(l) = tmp(l);
%   end
%   ocv(idx) = norm((f_r-f_e)./(1-h))^2;
%   gcv(idx) = norm((f_r-f_e)./(1-mean(h)))^2;
%  t_exact = t_exact+t_tmp+toc;
end
close(wb);


%% computations for plotting

[~,idx_err]       = min(err);
[~,idx_gcv]       = min(gcv);
[~,idx_ocv]       = min(ocv);
[~,idx_gcv_appr]  = min(gcv_appr);
[~,idx_ocv_appr]  = min(ocv_appr);

f_hat_r = compute_f_hat(plan_nfft,plan_infft,f_e,lambda(idx_ocv_appr),W_hat);

nfft_finalize(plan_nfft);

res = 480;
plotnodes = linspace(0,1,res);

plan_nfft = nfft_init_1d(N,res);
nfft_set_x(plan_nfft,plotnodes);
nfft_precompute_psi(plan_nfft);
plotf_r = F(plan_nfft,f_hat_r);
nfft_finalize(plan_nfft);


%% plotting

subplot(121);
% plot noisy data
scatter(nodes,real(f_e),5,'k','filled'); hold on;
% plot reconstruction
plot(plotnodes,real(plotf_r)); hold off;
axis square;
title('noisy data and reconstruction');

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

function f_hat = compute_f_hat(plan_nfft,plan_infft,f,lambda,W_hat)
  plan_infft.f = f;
  infft_trafo(plan_infft);
  f_hat = plan_infft.fcheck(:);
  [f_hat,~] = pcg(@(x) afun(plan_nfft,plan_infft,x,lambda,W_hat),f_hat);
end

function f_hat = compute_f_hat_appr(plan_nfft,plan_infft,y,lambda,W_hat)
  plan_infft.f = y;
  infft_trafo(plan_infft);
  f_hat = plan_infft.fcheck(:);
  f_hat = f_hat./(1+lambda*W_hat);
end

function f = F(plan,f_hat)
  nfft_set_f_hat(plan,f_hat);
  nfft_trafo(plan);
  f = nfft_get_f(plan);
end

function y = afun(plan_nfft,plan_infft,x,lambda,W_hat)
  nfft_set_f_hat(plan_nfft,x);
  nfft_trafo(plan_nfft);
  y = nfft_get_f(plan_nfft);

  plan_infft.f = y;
  infft_trafo(plan_infft);
  y = plan_infft.fcheck(:);
  
  y = y+lambda*W_hat.*x;
end




