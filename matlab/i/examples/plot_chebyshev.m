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
f      = f-min(f); f = f/max(f);      % normalize function

f_e    = f+0.05*randn(size(f));       % noisy function values
W      = 2/M;                         % weights in space domain      
s      = 3;

lambda = 2.^(linspace(-18,-11,25));   % possible lambda
err    = 0*lambda;                    % stores l_2-error
ocv    = 0*lambda;                    % stores ocv score
gcv    = 0*lambda;                    % stores gcv score


%% main computations

fcv = FCV_chebyshev(f_e,s);

for idx = 1:length(lambda) % loop over lambda
  [ocv(idx),gcv(idx),~,f_r] = fcv.compute(lambda(idx));
  err(idx) = norm(f-f_r);
end


%% computations for plotting

[~,idx_gcv] = min(gcv);
[~,idx_ocv] = min(ocv);
[~,~,~,f_r] = fcv.compute(lambda(idx_ocv));


%% plotting

subplot(121);
% plot noisy data
scatter(nodes,real(f_e),10,'k','filled'); hold on;
% plot reconstruction
plot(nodes,real(f_r)); hold off;
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
