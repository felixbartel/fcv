%
% this script calculates and plots the cv score for an equispaced example
% on the two-dimensional torus
%
%% initialization

rng('default');                                     % reset random generator

fun   = @(x,y) peaks(6*x-3,6*y-3);                  % example function
N     = 1024;                                       % bandwidth
[nodes_x,nodes_y] = meshgrid((0:N-1)/N,(0:N-1)/N);  % nodes in space domain
f     = reshape(fun(nodes_x,nodes_y),[],1);         % function values
f     = f-min(f); f = f/max(f);                     % normalize function
f_hat = 1/N^2*fft2(f,'transp');                     % get original f_hat

f_e   = f+0.1*randn(size(f));                       % noisy function values
W     = 1/N^2;                                      % weights in space domain
s     = 3;                                          % weights in frequency domain
W_hat = (1+[0:N/2-1 N/2:-1:1].^2'+[0:N/2-1 N/2:-1:1].^2).^(s/2);
W_hat = W_hat(:);

lambda = 2.^(linspace(-18,-8,25));                  % possible lambda
err    = 0*lambda;                                  % stores L_2-error
cv     = 0*lambda;                                  % stores cv score

t_slow = 0;
t_fast = 0;


%% main computations

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);

  tic;
  f_hat_r = compute_f_hat(f_e,lambda(idx),W,W_hat);
  f_r = F(f_hat_r);
  t_tmp = toc;
  err(idx) = norm(f_hat-f_hat_r);
  
% fast
  tic;
  h = W*sum(1./(1+lambda(idx)*W_hat));
  cv(idx) = norm((f_r-f_e)./(1-h))^2;
  t_fast = t_fast+t_tmp+toc;
  
% % slow
%   tic
%   h = zeros(N,1);
%   for l = 1:N^2
%     tmp = double( 1:N^2 == l );
%     tmp = F(compute_f_hat(tmp,lambda(idx),W,W_hat));
%     h(l) = tmp(l);
%   end
%   t_slow = t_slow+t_tmp+toc;
end
close(wb);

%% computations for plotting

[~,idx_err] = min(err);
[~,idx_cv]  = min(cv);

% calculate reconstruction
f_r = F(compute_f_hat(f_e,lambda(idx_cv),W,W_hat));

t_slow = t_slow/length(lambda);
t_fast = t_fast/length(lambda);


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


%% print times

fprintf('average time per lambda\n');
fprintf('\tt_slow = %f\n',t_slow);
fprintf('\tt_fast = %f\n',t_fast);


%% helper functions

function f_hat = compute_f_hat(y,lambda,W,W_hat)
  f_hat = W*fft2(y,'transp');
  f_hat = f_hat./(1+lambda*W_hat);
end

function f = F(f_hat)
  f = fft2(f_hat,'notransp');
end
