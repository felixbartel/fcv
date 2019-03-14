function [cv,grad,fhat_r,f_r] = compute_with_grad(self,lambda)
% fcv.COMPUTE_WITH_GRAD computes the cross-validation score
%
% Syntax:
%   cv = fcv.COMPUTE_WITH_GRAD(lambda)
%   [cv,grad,fhat_r,f_r] = fcv.COMPUTE_WITH_GRAD(lambda)
%
% Input:
%   lambda - regularization parameter 
%
% Output:
%   cv      - cross-validation score
%   grad    - gradient of the cross-validation score
%   fhat_r  - corresponding Fourier coefficients
%   f_r     - corresponding function values

  fhat_r = self.W*fftd_adj(self.f,self.d);
  fhat_r = fhat_r./(1+lambda*self.What);
  f_r = fftd(fhat_r,self.d);

  h = self.W*sum(1./(1+lambda*self.What));
  cv = sum(((f_r-self.f)./(1-h)).^2);
  cv = real(cv);
  
  ddh = self.W*sum(self.What./(1+lambda*self.What).^2);
  ddf = fftd(fhat_r.*self.What./(1+lambda*self.What).^2,self.d);
  
  grad = -2*sum(...
    ddh.*(f_r-self.f).^2./(1-h).^3+...
    ddf.*(f_r-self.f)./(1-h).^2);
  grad = real(grad);
end