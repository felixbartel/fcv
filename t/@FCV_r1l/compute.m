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

  fhat_r = alfft(self.z,self.M,self.I,self.W*self.f);
  fhat_r = fhat_r./(1+lambda*self.What);
  
  f_r = lfft(self.z,self.M,self.I,fhat_r.');
  h = sum(1./(1+lambda*self.What))*self.W;
  cv = 1/self.M*real(sum(((f_r-self.f)./(1-h)).^2));  
end
