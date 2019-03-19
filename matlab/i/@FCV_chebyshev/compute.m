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
  fhat_r(1) = sqrt(2)*fhat_r(1);
  fhat_r = sqrt(self.N/2)*fhat_r;
  
  fhat_r = fhat_r./(pi*[1; 1/2*ones(self.N-1,1)]+lambda*self.What);
  
  fhat_r(1) = sqrt(2)*fhat_r(1);
  f_r = sqrt(self.N/2)*dctIII(fhat_r);

% compute diagonal emelents via dctI
  b = 1./(pi*[1; 1/2*ones(self.N-1,1)]+lambda*self.What);

  btilde = zeros(2*self.M+1,1);
  btilde(1:2:2*self.M-1) = b;
  btilde(1) = btilde(1)*2*sqrt(2);

  h = sqrt(self.N/2)*dctI(btilde);
  h = h(2:2:2*self.M);
  h = h+(sum(b(2:end)));
  
  h = self.W/2*h;
  
  ocv = norm((f_r-self.f)./(1-h))^2;
  gcv = norm((f_r-self.f)./(1-mean(h)))^2;
end