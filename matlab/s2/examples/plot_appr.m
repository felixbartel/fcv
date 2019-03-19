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
ocv_exact = 0*lambda;                     % stores ocv score
gcv_exact = 0*lambda;                     % stores gcv score
ocv       = 0*lambda;                     % stores approximated ocv score
gcv       = 0*lambda;                     % stores approximated gcv score


%% original fhat

sF = S2FunHarmonic.quadrature(@(v) fun(v));
fhat = sF.fhat(1:(N+1)^2);


%% main computations

fcv = FCV_appr([nodes.rho,nodes.theta],f_e,W,N,s);

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
  [ocv(idx),gcv(idx),fhat_r] = fcv.compute(lambda(idx));
  err(idx) = norm(fhat-fhat_r);
%  [ocv_exact(idx),gcv_exact(idx),~] = fcv.compute_exact(lambda(idx));
end
close(wb);

%% computations for plotting

[~,idx_err]       = min(err);
[~,idx_gcv_exact] = min(gcv_exact);
[~,idx_ocv_exact] = min(ocv_exact);
[~,idx_gcv]       = min(gcv);
[~,idx_ocv]       = min(ocv);

[~,~,fhat_r] = fcv.compute(lambda(idx_ocv));
sF_r = S2FunHarmonic(fhat_r);


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
p = loglog(lambda,[ocv_exact; gcv_exact; ocv; gcv]); hold on;
scatter(lambda([idx_ocv_exact idx_gcv_exact idx_ocv idx_gcv]),...
  [ocv_exact(idx_ocv_exact) gcv_exact(idx_gcv_exact) ocv(idx_ocv) gcv(idx_gcv)],40,'filled');
ylim([min([ocv_exact gcv_exact gcv]) max([ocv_exact gcv_exact gcv])]); hold off;
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