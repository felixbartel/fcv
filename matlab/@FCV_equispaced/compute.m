function [cv,fhat_r,f_r] = compute(self,lambda)
% fcv.COMPUTE computes the cross-validation score
%
% Syntax:
%   cv = fcv.COMPUTE(lambda)
%   [cv,f_hat_r,f_r] = fcv.COMPUTE(lambda)
%
% Input:
%   lambda - regularization parameter 
%
% Output:
%   cv      - cross-validation score
%   f_hat_r - corresponding Fourier coefficients
%   f_r     - corresponding function values

  fhat_r = self.W*fftd_adj(self.f,self.d);
  fhat_r = fhat_r./(1+lambda*self.What);
  f_r = fftd(fhat_r,self.d);

  h = self.W*sum(1./(1+lambda*self.What));
  cv = sum(((f_r-self.f)./(1-h)).^2);
  cv = real(cv);
end