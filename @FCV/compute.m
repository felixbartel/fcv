function s = compute(self,lambda,varargin)
% fcv.COMPUTE computes the approximated cross-validation score
%
% Syntax:
%   s = fcv.COMPUTE(lambda)
%   s = fcv.COMPUTE(lambda,'exact')
%   s = fcv.COMPUTE(lambda,'derivative')
%
% Input:
%   lambda - regularization parameter 
%
% Output:
%   s.ocv     - ordinary cross-validation score
%   s.gcv     - generalized cross-validation score
%   s.fhat_r  - corresponding Fourier coefficients
%   s.f_r     - corresponding function values

  flag_exact = "";
  derivative = false;
  for j = varargin
    if strcmp(j{:},"exact"); flag_exact = "exact"; end
    if strcmp(j{:},"derivative"); derivative = true; end
  end
  if derivative
    [s.f_r,s.fhat_r,ddf_r] = self.H(lambda);
    [h,ddh] = self.diagonals(lambda,flag_exact);
    s.ddocv = 2/self.M*sum(...
      ddh.*(s.f_r-self.f).^2./(1-h).^3+...
      ddf_r.*(s.f_r-self.f)./(1-h).^2);
    s.ddgcv = 2/self.M*sum(...
      mean(ddh).*(s.f_r-self.f).^2./(1-mean(h)).^3+...
      ddf_r.*(s.f_r-self.f)./(1-mean(h)).^2);
  else
    [s.f_r,s.fhat_r] = self.H(lambda);
    h = self.diagonals(lambda,flag_exact);
  end
  s.ocv = 1/self.M*norm((s.f_r-self.f)./(1-h)).^2;
  s.gcv = 1/self.M*norm((s.f_r-self.f)./(1-mean(h))).^2;
end
