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
s         = 3;

lambda    = 2.^linspace(-17,-11,25);        % possible lambda
err       = 0*lambda;                       % stores l_2-error
ocv_exact = 0*lambda;                       % stores ocv score
gcv_exact = 0*lambda;                       % stores gcv score
ocv       = 0*lambda;                       % stores approximated ocv score
gcv       = 0*lambda;                       % stores approximated gcv score


%% main computations

fcv = FCV_appr(nodes,f_e,[],M,s);

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
  [ocv(idx),gcv(idx),~,f_r] = fcv.compute(lambda(idx));
  err(idx) = norm(f-f_r);
%  [ocv_exact(idx),gcv_exact(idx),~,f_r] = fcv.compute_exact(lambda(idx));
end
close(wb);


%% computations for plotting

[~,idx_gcv_exact] = min(gcv_exact);
[~,idx_ocv_exact] = min(ocv_exact);
[~,idx_gcv] = min(gcv);
[~,idx_ocv] = min(ocv);

[~,~,f_hat_r] = fcv.compute(lambda(idx_ocv));

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
p = loglog(lambda,[ocv_exact;gcv_exact;ocv;gcv]); hold on;
scatter(lambda([idx_ocv_exact idx_gcv_exact idx_ocv idx_gcv]),...
  [ocv_exact(idx_ocv_exact) gcv_exact(idx_gcv_exact) ocv(idx_ocv) gcv(idx_gcv)],40,'filled'); hold off;
ylim([min([ocv_exact gcv_exact gcv]) max([ocv_exact gcv_exact gcv])]); hold off;
legend(p,'ocv','gcv','appr ocv','appr gcv');
ylabel('cv score');
axis square;
xlim([lambda(1) lambda(end)]);
