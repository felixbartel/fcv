function [ocv,gcv,fhat_r,f_r] = compute_exact(self,lambda)
% fcv.COMPUTE_EXACT computes the exact cross-validation score
%
% Syntax:
%   ocv = fcv.COMPUTE_EXACT(lambda)
%   [ocv,gcv,f_hat_r,f_r] = fcv.COMPUTE_EXACT(lambda)
%
% Input:
%   lambda - regularization parameter 
%
% Output:
%   ocv     - ordinary cross-validation score
%   gcv     - generalized cross-validation score
%   f_hat_r - corresponding Fourier coefficients
%   f_r     - corresponding function values

  [fhat_r,~] = lsqr(...
    @(x,transp_flag) A(self.plan,x,lambda,self.W,self.What,transp_flag),...
    [sqrt(self.W).*self.f;zeros(length(self.What),1)]);
    
  nfft_set_f_hat(self.plan,fhat_r);
  nfft_trafo(self.plan);
  f_r = nfft_get_f(self.plan);

% exact cv score
  h = zeros(self.M,1);
  for l = 1:self.M
    tmp = double( 1:self.M == l )';
    
    [fhat_tmp,~] = lsqr(...
    @(x,transp_flag) A(self.plan,x,lambda,self.W,self.What,transp_flag),...
    [sqrt(self.W).*tmp;zeros(length(self.What),1)]);
    
    nfft_set_f_hat(self.plan,fhat_tmp);
    nfft_trafo(self.plan);
    tmp = nfft_get_f(self.plan);
    
    h(l) = tmp(l);
  end
  ocv = real(sum(((f_r-self.f)./(1-h)).^2));  
  gcv = real(sum(((f_r-self.f)./(1-mean(h))).^2));
end


function y = A(plan,x,lambda,W,What,transp_flag)
  if strcmp(transp_flag,'notransp')
    nfft_set_f_hat(plan,x);
    nfft_trafo(plan);
    y = nfft_get_f(plan);
    
    y = sqrt(W).*y;
    
    y = [y; sqrt(lambda*What).*x];
  else
    y = sqrt(W).*x(1:length(W));
    
    nfft_set_f(plan,y);
    nfft_adjoint(plan);
    y = nfft_get_f_hat(plan);
    
    y = y+sqrt(lambda*What).*x(length(W)+1:end);
  end
end
