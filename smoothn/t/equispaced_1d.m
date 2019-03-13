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


lambda_0 = 1;
MaxFunEvals = 10;

%% plot noisy data

subplot(121);
scatter(nodes,real(f_e),5,'k','filled'); hold on;
axis square;
title('noisy data');


%% main computations

options = optimoptions('fminunc',...
  'SpecifyObjectiveGradient',true,...
  'OutputFcn',@outfun,...
  'Display','iter',...
  'MaxFunEvals',MaxFunEvals);

subplot(122);
lambda_min = fminunc(@(t) cvfun(f_e, t, W, W_hat),...
  log(lambda_0), options);
hold off;
lambda_min = exp(lambda_min);



%% calculate reconstruction
f_r = F(compute_f_hat(f_e,lambda_min,W,W_hat));


%% plot reconstruction

subplot(121);
plot(nodes,real(f_r)); hold off;
title('noisy data and reconstruction');


%% helper functions

function stop = outfun(x, optimValues, ~)
  stop = false;
  loglog(exp(x),exp(optimValues.fval),'k.');
  xlabel('\lambda');
  ylabel('cv score');
  axis square;
  hold on;
  drawnow;
end

function [val,grad] = cvfun(f,t,W,W_hat)
  f_hat = W*length(f)*ifft(f);
  f_1 = F(f_hat./(1+exp(t)*W_hat));
  f_2 = F(f_hat.*W_hat./(1+exp(t)*W_hat).^2);
  
  h = W*sum(1./(1+exp(t)*W_hat));
  val = sum(((f_1-f)./(1-h)).^2);
  val = real(val);
  
  h2 = W*sum(W_hat./(1+exp(t)*W_hat).^2);
  
  grad = -2*sum(...
    h2.*(f_1-f).^2./(1-h).^3+...
    f_2.*(f_1-f)./(1-h).^2);
  grad = real(grad)/val*exp(t);
  val = log(val);
end

% variant not in log scale 
% function [val,grad] = cvfun(f,lambda,W,W_hat)
%   f_hat = W*length(f)*ifft(f);
%   f_1 = F(f_hat./(1+lambda*W_hat));
%   f_2 = F(f_hat.*W_hat./(1+lambda*W_hat).^2);
%   
%   h = W*sum(1./(1+lambda*W_hat));
%   val = sum(((f_1-f)./(1-h)).^2);
%   val = real(val);
%   
%   h2 = W*sum(W_hat./(1+lambda*W_hat).^2);
%   
%   grad = -2*sum(...
%     h2.*(f_1-f).^2./(1-h).^3+...
%     f_2.*(f_1-f)./(1-h).^2);
%   grad = real(grad);
% end

function f_hat = compute_f_hat(f,lambda,W,W_hat)
  f_hat = W*length(f)*ifft(f);
  f_hat = f_hat./(1+lambda*W_hat);
end

function f = F(f_hat)
  f = fft(f_hat);
end