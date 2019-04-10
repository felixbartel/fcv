K%
% this script calculates and plots the cv score for an equispaced example
% on the one-dimensional torus
%
%% initialization

rng('default');                      % reset random generator

fun    = @(x) peaks(6*x-3,0);        % example function
M      = 2^10;                       % number of nodes
nodes  = (0:M-1)'/M;                 % nodes in space domain
f      = fun(nodes);                 % function values
f = f-min(f); f = f/max(f);          % normalize function
fhat   = ifft(f);                    % get original fhat

f_e    = f+0.1*randn(size(f));       % noisy function values
s      = 3;                          % weights in frequency domain

lambda = 2.^(linspace(-20,-5,25));   % possible lambda
err    = 0*lambda;                   % stores L_2-error
cv     = 0*lambda;                   % stores cv score


%% main computations

fcv = FCV_equispaced(1,f_e,s);

for idx = 1:length(lambda) % loop over lambda
  [cv(idx),f_hat_r] = fcv.compute(lambda(idx));
  err(idx) = norm(fhat-f_hat_r);
end


%% computations for plotting

[~,idx_err] = min(err);
[~,idx_cv]  = min(cv);

% calculate reconstruction
[~,~,f_r] = fcv.compute(lambda(idx_cv));


%% plotting

subplot(121);
% plot noisy data
scatter(nodes,real(f_e),5,'k','filled'); hold on;
% plot reconstruction
plot(nodes,real(f_r)); hold off;
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
loglog(lambda,cv); hold on;
scatter(lambda(idx_cv),cv(idx_cv),'filled'); hold off;
ylabel('cv score');
axis square;




