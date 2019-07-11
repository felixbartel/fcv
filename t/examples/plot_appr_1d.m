%
% this script calculates and plots the cv score for an nonequispaced
% example on the one-dimensional torus
%
%% initialization

rng('default');                       % reset random generator
fun       = @(x) peaks(6*x-3,0);      % example function 
N         = 2^6;                      % bandwidth
nodes     = rand(2*N,1);              % nodes in space domain
nodes     = nodes.^2;
nodes     = unique(nodes,'rows');
f         = fun(nodes);               % function values
f = f-min(f); f = f/max(f);           % normalize function

f_e       = f+0.05*randn(size(f));    % noisy function values
s         = 3;                        % weights in frequency domain

lambda    = 2.^(linspace(-15,-5,25)); % possible lambda
err       = 0*lambda;                 % stores L_2-error
ocv_exact = 0*lambda;                 % stores ocv score
gcv_exact = 0*lambda;                 % stores gcv score
ocv       = 0*lambda;                 % stores approximated ocv score
gcv       = 0*lambda;                 % stores approximated gcv score


%% original f_hat

M2 = 2^10;                          % sufficiently large bandwidth
[nodes2] = (0:M2-1)/M2;             % nodes in space domain
fe = fun(nodes2);                   % function values
fe = fe-min(fe); fe = fe/max(fe);
fhat = ifft(fe);
fhat = fhat([(end-N/2+1):end 1:N/2]);


%% main computations

fcv = FCV_appr(nodes,f_e,[],N,s);

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
  s = fcv.compute(lambda(idx));
  ocv(idx) = s.ocv;
  gcv(idx) = s.gcv;
%   s = fcv.compute(lambda(idx),"exact");
%   ocv_exact(idx) = s.ocv;
%   gcv_exact(idx) = s.gcv;
  err(idx) = norm(fhat-s.fhat_r);
end
close(wb);


%% computations for plotting

[~,idx_err]       = min(err);
[~,idx_gcv]       = min(gcv_exact);
[~,idx_ocv]       = min(ocv_exact);
[~,idx_gcv_appr]  = min(gcv);
[~,idx_ocv_appr]  = min(ocv);

s = fcv.compute(lambda(idx_ocv_appr));

res = 480;
plotnodes = linspace(0,1,res);

plan = nfft_init_1d(N,res);
nfft_set_x(plan,plotnodes);
nfft_precompute_psi(plan);
nfft_set_f_hat(plan,s.fhat_r);
nfft_trafo(plan);
plotf_r = nfft_get_f(plan);
nfft_finalize(plan);


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
p = loglog(lambda,[ocv_exact; gcv_exact; ocv; gcv]); hold on;
scatter(lambda([idx_ocv idx_gcv idx_ocv_appr idx_gcv_appr]),...
  [ocv_exact(idx_ocv) gcv_exact(idx_gcv) ocv(idx_ocv_appr) gcv(idx_gcv_appr)],40,'filled');
ylim([min(real([ocv_exact gcv_exact gcv])) max(real([ocv_exact gcv_exact gcv]))]); hold off;
legend(p,'ocv','gcv','appr ocv','appr gcv');
ylabel('cv score');
axis square;






