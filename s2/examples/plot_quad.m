%
% this script calculates and plots the cv score for an quadrature
% example on the two-dimensional sphere
%
%% initialization

rng('default');                           % reset random generatormber of nodes
fun       = @(nodes) S2Fun.smiley(nodes); % example function
N         = 100;                          % bandwidth
[nodes,W] = quadratureS2Grid(2*N);        % nodes and weights in space domain
M         = length(nodes);                % number of nodes
f         = fun(nodes);                   % function values

f_e       = f+0.1*randn(size(f));         % noisy function values
s         = 3;                            % weights in requency domain

lambda    = 2.^(linspace(-38,-25,25));    % possible lambda
err       = 0*lambda;                     % stores L_2-error
cv        = 0*lambda;                     % stores cv score


%% original fhat

sF = S2FunHarmonic.quadrature(@(v) fun(v));
fhat = sF.fhat(1:(N+1)^2);


%% main computations

fcv = FCV_quad([nodes.rho,nodes.theta],f_e,W,N,s);

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
  res = fcv.compute(lambda(idx));
  cv(idx) = res.ocv;
  err(idx) = norm(fhat-res.fhat_r);
end
close(wb);


%% computations for plotting

% get optimal lambda
[~,idx_err] = min(err);
[~,idx_cv]  = min(cv);

res = fcv.compute(lambda(idx_cv));
sF_r = S2FunHarmonic(res.fhat_r);


%% plotting

figure(1);
% plot L_2-error
yyaxis left;
loglog(lambda,err); hold on;
scatter(lambda(idx_err),err(idx_err),40,'filled');
xlabel('\lambda');
ylabel('L_2-error');
% plot cv score
yyaxis right;
p = loglog(lambda,cv); hold on;
scatter(lambda(idx_cv),cv(idx_cv),40,'filled');
ylabel('cv score');
axis square;

% plot functions on the sphere
figure(2);
% plot noisy data
pcolor(nodes,f_e,'contours',20,'upper');
setColorRange([-0.5 0.5]);
%hold on; scatter(nodes,'MarkerSize',2,'MarkerColor','k'); hold off
mtexTitle('noisy data');
nextAxis;
% plot reconstruction
pcolor(sF_r,'contours',20,'upper');
setColorRange([-0.5 0.5]);
mtexTitle('reconstruction');
mtexColorbar;
