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

  self.plan.f = self.W.*self.f;
  nfsoft_adjoint(self.plan);
  fhat_r = self.plan.fhat;

  tmp = 2*(0:self.N)+1;
  tmp = repelem(tmp,(1:2:(2*self.N+1)).^2)';
  
  fhat_r = fhat_r./(8*pi^2./tmp+lambda*self.What);
  self.plan.fhat = fhat_r;
  nfsoft_trafo(self.plan);
  f_r = self.plan.f;
  
  idx = nfsoft_f_hat_size(0:self.N);
  h = sum(tmp(idx)./(8*pi^2./tmp(idx)+lambda*self.What(idx)))*self.W;
  ocv = real(sum(((f_r-self.f)./(1-h)).^2));
  gcv = real(sum(((f_r-self.f)./(1-mean(h))).^2));
end
