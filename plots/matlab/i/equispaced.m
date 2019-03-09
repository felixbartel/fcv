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
nu     = 3;
W_hat  = ((1:M).').^nu;               % weights in frequency domain

lambda = 2.^(linspace(-18,-11,25));   % possible lambda
err    = 0*lambda;                    % stores l_2-error
ocv    = 0*lambda;                    % stores ocv score
gcv    = 0*lambda;                    % stores gcv score


%% main computations

for idx = 1:length(lambda) % loop over lambda
  f_r = mult_H(f_e,lambda(idx),W,W_hat);
  err(idx) = norm(f-f_r);

% compute diagonal emelents via dctI
  a = 1./(1+lambda(idx)*W_hat);

  atilde = zeros(2*M+1,1);
  atilde(1:2:2*M-1) = a;
  atilde(1) = atilde(1)/sqrt(2);

  h = sqrt(1/M)*dctI(atilde);
  h = h(2:2:2*M);
  h = h+(sum(a)-a(1)/2)/M;

  ocv(idx) = norm((f_r-f_e)./(1-h))^2;
  gcv(idx) = norm((f_r-f_e)./(1-mean(h)))^2;
end


%% computations for plotting

[~,idx_gcv] = min(gcv);
[~,idx_ocv] = min(ocv);
f_r = mult_H(f_e,lambda(idx_ocv),W,W_hat);


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


%% helper functions

function f = mult_H(f,lambda,W,W_hat)
  M = length(f);
%F = sqrt(2/N)*cos((pi*(2*(0:N-1)+1)/(2*N)).'*(0:N-1)); % dctIII
%F(:,1) = F(:,1)/sqrt(2);

% multiplies with the hat matrix
  f_hat = dctII(W*f);
  f_hat = f_hat./(1+lambda*W_hat);
  f = M/2*dctIII(f_hat);
end

function ahat = dctI(a)
  s = size(a);
  a = a(:);
  N = length(a);

  y = a([1:N N-1:-1:2]);
  y([1 N]) = sqrt(2)*y([1 N]);

  yhat = fft(y);

  ahat = 1/sqrt(2*N-2)*real(yhat(1:N));
  ahat([1 end]) = ahat([1 end])/sqrt(2);

  reshape(ahat,s);
end

function ahat = dctII(a)
  s = size(a);
  a = a(:);
  N = length(a);

  y = a([1:N N:-1:1]);

  yhat = fft(y);

  ahat = 1/sqrt(2*N)*real(exp(-2i*pi*(0:N-1)'/(4*N)).*yhat(1:N));
  ahat(1) = ahat(1)/sqrt(2);

  reshape(ahat,s);
end

function ahat = dctIII(a)
  s = size(a);
  a = a(:);
  N = length(a);

  y = [sqrt(2)*a(1); ...
    exp(-2i*pi*(1:N-1)'/(4*N)).*a(2:N); ...
    0; ...
    exp(-2i*pi*((3*N+1):(4*N-1))'/(4*N)).*a(N:-1:2)];
  yhat = fft(y);
  ahat = 1/sqrt(2*N)*real(yhat(1:N));

  ahat = reshape(ahat,s);
end
