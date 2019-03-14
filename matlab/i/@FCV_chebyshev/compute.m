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

  fhat_r = dctII(self.W*self.f);
  fhat_r = fhat_r./(1+lambda*self.What);
  f_r = self.M/2*dctIII(fhat_r);

% compute diagonal emelents via dctI
  a = 1./(1+lambda*self.What);

  atilde = zeros(2*self.M+1,1);
  atilde(1:2:2*self.M-1) = a;
  atilde(1) = atilde(1)/sqrt(2);

  h = sqrt(1/self.M)*dctI(atilde);
  h = h(2:2:2*self.M);
  h = h+(sum(a)-a(1)/2)/self.M;

  ocv = norm((f_r-self.f)./(1-h))^2;
  gcv = norm((f_r-self.f)./(1-mean(h)))^2;
end