%
% this script calculates and plots the cv score for an equispaced example
% on the two-dimensional torus
%
%% initialization

rng('default');                                     % reset random generator

fun    = @(x,y) peaks(6*x-3,6*y-3);                 % example function
N      = 2^10;                                      % bandwidth
[nodes_x,nodes_y] = meshgrid((0:N-1)/N,(0:N-1)/N);  % nodes in space domain
f      = reshape(fun(nodes_x,nodes_y),[],1);        % function values
f = f-min(f); f = f/max(f);                         % normalize function
fhat   = fftd_adj(f,2)/N^2;                         % get original fhat

f_e    = f+0.1*randn(size(f));                      % noisy function values
s      = 3;                                         % weights in frequency domain

lambda = 2.^(linspace(-18,-8,25));                  % possible lambda
err    = 0*lambda;                                  % stores L_2-error
cv     = 0*lambda;                                  % stores cv score


%% main computations

fcv = FCV_equispaced(2,f_e,s);

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
  [cv(idx),f_hat_r] = fcv.compute(lambda(idx));
  err(idx) = norm(fhat-f_hat_r);
end
close(wb);

%% computations for plotting

[~,idx_err] = min(err);
[~,idx_cv]  = min(cv);

% calculate reconstruction
[~,~,f_r] = fcv.compute(lambda(idx_cv));


%% plotting

% plot noisy data
subplot(221);
imagesc(real(reshape(f_e,N,N)));
title('noisy data');
axis square;
caxis([0 1]);
colorbar;

% plot reconstruction
subplot(223);
imagesc(real(reshape(f_r,N,N)));
title('reconstruction');
axis square;
caxis([0 1]);
colorbar;

colormap jet;

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
