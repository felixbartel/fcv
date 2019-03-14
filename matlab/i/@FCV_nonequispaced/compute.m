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

  f = sqrt(self.M/2*self.W).*self.f;
  f = [f;zeros(self.M,1)];
  [fhat_r,~] = lsqr(...
    @(x,transp_flag) afun(self.plan,x,lambda,self.W,self.What,transp_flag),...
    f,1e-8,1000);
  
  f_r = ndctIII(self.plan,fhat_r); 
  
% compute diagonal emelents
  a = 1./(1+lambda*self.What);
  a(1) = a(1)/2;

  a_hat_r = zeros(2*self.M,1);
  a_hat_r(1:2:2*self.M-1) = a;

  nfct_set_f_hat(self.plan2,double(a_hat_r));
  nfct_trafo(self.plan2);
  h = nfct_get_f(self.plan2);

  h = self.W.*(h+sum(a))/2;

  ocv = real(sum(((f_r-self.f)./(1-h)).^2));
  gcv = real(sum(((f_r-self.f)./(1-mean(h))).^2));
end


function y = afun(plan,x,lambda,W,What,transp_flag)
  M = length(W);
  if strcmp(transp_flag,'notransp')
    y = sqrt(M/2*W).*ndctIII(plan,x);
    y = [y; sqrt(lambda*What).*x];
  else
    y = ndctII(plan,sqrt(M/2*W).*x(1:M));
    y = y+sqrt(lambda*What).*x(M+1:end);
  end
end


