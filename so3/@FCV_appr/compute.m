function [ocv,gcv,fhat_r,f_r] = compute(self,lambda)
% fcv.COMPUTE computes the cross-validation score
%
% Syntax:
%   ocv = fcv.COMPUTE(lambda)
%   [ocv,gcv,fhat_r,f_r] = fcv.COMPUTE(lambda)
%
% Input:
%   lambda - regularization parameter 
%
% Output:
%   ocv     - ordinary cross-validation score
%   gcv     - generalized cross-validation score
%   fhat_r  - corresponding Fourier coefficients
%   f_r     - corresponding function values

  [fhat_r,~] = lsqr(...
  @(x,transp_flag) A(self.plan,x,lambda,self.W,self.What,transp_flag),...
  [sqrt(self.W).*self.f;zeros(length(self.What),1)],1e-10);
  
  self.plan.fhat = fhat_r;
  nfsoft_trafo(self.plan);
  f_r = self.plan.f;
  
  tmp = 2*(0:self.N)+1;
  tmp = repelem(tmp,(1:2:(2*self.N+1)).^2)';
  
  idx = nfsoft_f_hat_size(0:self.N);
  h = sum(tmp(idx)./(8*pi^2./tmp(idx)+lambda*self.What(idx)))*self.W;
  ocv = 1/self.M*real(sum(((f_r-self.f)./(1-h)).^2));
  gcv = 1/self.M*real(sum(((f_r-self.f)./(1-mean(h))).^2));
end



function y = A(plan,x,lambda,W,W_hat,transp_flag)
  if strcmp(transp_flag,'notransp')
    plan.fhat = x;
    nfsoft_trafo(plan);
    y = plan.f;
    
    y = sqrt(W).*y;
    
    y = [y; sqrt(lambda*W_hat).*x];
  else
    y = sqrt(W).*x(1:length(W));
    
    plan.f = y;
    nfsoft_adjoint(plan);
    y = plan.fhat;
    
    y = y+sqrt(lambda*W_hat).*x(length(W)+1:end);
  end
end
