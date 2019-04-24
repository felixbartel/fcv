%
% this script calculates and plots the cv score for an nonequispaced
% example on the two-dimensional torus
%
%% initialization

rng('default');                        % reset random generator
fun       = @(x,y) peaks(6*x-3,6*y-3); % example function 
N         = 64;                        % bandwidth
nodes     = rand(2*N^2,2);             % nodes in space domain
nodes     = nodes.^2;
nodes     = uniquetol(nodes,1e-5,'ByRows',true);
f         = fun(nodes(:,1),nodes(:,2));% function values
f = f-min(f); f = f/max(f);

f_e       = f+0.05*randn(size(f));     % noisy function values
s         = 3;                         % weights in frequency domain

lambda    = 2.^(linspace(-12,-6,25));  % possible lambda
err       = 0*lambda;                  % stores L_2-error
ocv_exact = 0*lambda;                  % stores ocv score
gcv_exact = 0*lambda;                  % stores approximated ocv score
ocv       = 0*lambda;                  % stores gcv score
gcv       = 0*lambda;                  % stores approximated gcv score


%% original f_hat

M2 = 2^10;                                              % sufficiently large bandwidth
[nodes_x,nodes_y] = meshgrid((0:M2-1)/M2,(0:M2-1)/M2);  % nodes in space domain
fe = reshape(fun(nodes_x,nodes_y),[],1);                % function values
fe = fe-min(fe); fe = fe/max(fe);

fhat = 1/M2^2*fftd_adj(fe,2);
fhat = reshape(fhat,M2,M2);
fhat = fhat([(M2-N/2+1):M2 1:(N/2)],[(M2-N/2+1):M2 1:(N/2)]);
fhat = fhat(:);


%% main computations

fcv = FCV_appr(nodes,f_e,[],N,s);

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
  [ocv(idx),gcv(idx),f_hat_r] = fcv.compute(lambda(idx));
%  [ocv_exact(idx),gcv_exact(idx)] = fcv.compute_exact(lambda(idx));
  err(idx) = norm(fhat-f_hat_r);
end
close(wb);


%% computations for plotting

[~,idx_err]       = min(err);
[~,idx_gcv]       = min(gcv_exact);
[~,idx_ocv]       = min(ocv_exact);
[~,idx_gcv_appr]  = min(gcv);
[~,idx_ocv_appr]  = min(ocv);

[~,~,f_hat_r] = fcv.compute(lambda(idx_ocv_appr));

res = 480;
t = linspace(0,1,res);
[plotnodes_x,plotnodes_y] = meshgrid(t,t);
plan = nfft_init_2d(N,N,res^2);
nfft_set_x(plan,[plotnodes_x(:).';plotnodes_y(:).']);
nfft_precompute_psi(plan);
nfft_set_f_hat(plan,f_hat_r);
nfft_trafo(plan);
plotf_r = nfft_get_f(plan);
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
p = loglog(lambda,[ocv_exact; gcv_exact; ocv; gcv]); hold on;
scatter(lambda([idx_ocv idx_gcv idx_ocv_appr idx_gcv_appr]),...
  [ocv_exact(idx_ocv) gcv_exact(idx_gcv) ocv(idx_ocv_appr) gcv(idx_gcv_appr)],40,'filled');
ylim([min([ocv_exact gcv_exact gcv]) max([ocv_exact gcv_exact gcv])]); hold off;
legend(p,'ocv','gcv','appr ocv','appr gcv');
ylabel('cv score');
axis square;
