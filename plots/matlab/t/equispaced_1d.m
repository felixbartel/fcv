%
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
f_hat = ifft(f);                     % get original f_hat

f_e    = f+0.1*randn(size(f));       % noisy function values
W      = 1/M;                        % weights in space domain
s      = 3;                          % weights in frequency domain
W_hat  = ([0:M/2-1 M/2:-1:1])'.^s+1;

lambda = 2.^(linspace(-20,-5,25));   % possible lambda
err    = 0*lambda;                   % stores L_2-error
cv     = 0*lambda;                   % stores cv score


%% main computations

for idx = 1:length(lambda) % loop over lambda
  f_hat_r = compute_f_hat(f_e,lambda(idx),W,W_hat);
  err(idx) = norm(f_hat-f_hat_r);
  f_r = F(f_hat_r);
  
  h = W*sum(1./(1+lambda(idx)*W_hat));
  cv(idx) = norm((f_r-f_e)./(1-h))^2;
end


%% computations for plotting

[~,idx_err] = min(err);
[~,idx_cv]  = min(cv);

% calculate reconstruction
f_r = F(compute_f_hat(f_e,lambda(idx_cv),W,W_hat));


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


%% helper functions

function f_hat = compute_f_hat(f,lambda,W,W_hat)
  f_hat = W*length(f)*ifft(f);
  f_hat = f_hat./(1+lambda*W_hat);
end

function f = F(f_hat)
  f = fft(f_hat);
end