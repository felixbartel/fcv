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

  nfsftmex('set_f',self.plan,self.W.*self.f);
  nfsftmex('adjoint',self.plan);
  fhat_r = nfsftmex('get_f_hat_linear',self.plan);

  fhat_r = fhat_r./(1+lambda*self.What);
  nfsftmex('set_f_hat_linear',self.plan,fhat_r);
  nfsftmex('trafo',self.plan);
  f_r = nfsftmex('get_f',self.plan);
  
  h = sum((2*(0:self.N)'+1)./(1+lambda*self.What((1:self.N+1).^2)))/(4*pi)*self.W;
  ocv = real(sum(((f_r-self.f)./(1-h)).^2));
  gcv = real(sum(((f_r-self.f)./(1-mean(h))).^2));
end
