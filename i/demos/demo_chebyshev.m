%
% this script calculates and plots the cv score for an example on the unit
% interval with chebyshev nodes
%
%% initialization

rng('default');                       % reset random generator
fun    = @(x) peaks(2*x,0);           % example function
M      = 2^7;                         % number of nodes
nodes  = cos(pi*(2*(1:M)-1)/(2*M)).'; % Tschebyshev nodes of first kind
f      = fun(nodes);                  % function values
f = f-min(f); f = f/max(f);           % normalize function

f_e    = f+0.05*randn(size(f));       % noisy function values
s      = 3;

lambda = 2.^(linspace(-18,-10,20));   % possible lambda
err    = 0*lambda;                    % stores l_2-error
ocv    = 0*lambda;                    % stores ocv score
gcv    = 0*lambda;                    % stores gcv score


%% compute exact fhat

M2 = 2^10;
nodes2 = cos(pi*(2*(1:M2)-1)/(2*M2)).';
f2 = fun(nodes2);
f2 = f2-min(f2); f2 = f2/max(f2);
fcv = FCV_chebyshev(f2,0);
res = fcv.compute(0);
fhat = res.fhat_r;
fhat = fhat(1:M);


%% main computations

fcv = FCV_chebyshev(f_e,s);

for idx = 1:length(lambda) % loop over lambda
  res = fcv.compute(lambda(idx));
  ocv(idx) = res.ocv;
  gcv(idx) = res.gcv;
  err(idx) = norm([pi;pi/2*ones(M-1,1)].*(fhat-res.fhat_r));
end


%% computations for plotting

[~,idx_gcv] = min(gcv);
[~,idx_ocv] = min(ocv);
res = fcv.compute(lambda(idx_ocv));


%% plotting

subplot(121);
% plot noisy data
scatter(nodes,real(f_e),10,'k','filled'); hold on;
% plot reconstruction
plot(nodes,real(res.f_r)); hold off;
axis square;
title('noisy data and recontruction')

subplot(122);
% plot l_2-error
yyaxis left;
loglog(lambda,err);
xlabel('\lambda');
ylabel('l_2-error');
% plot cv score
yyaxis right;
p = loglog(lambda,[ocv; gcv]); hold on;
scatter(lambda([idx_ocv idx_gcv]),...
  [ocv(idx_ocv) gcv(idx_gcv)],40,'filled'); hold off;
legend(p,'ocv','gcv');
ylabel('cv score');
axis square;
