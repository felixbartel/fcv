function [ocv,gcv,fhat_r,f_r] = compute_exact(self,lambda)
% fcv.COMPUTE_EXACT computes the exact cross-validation score
%
% Syntax:
%   ocv = fcv.COMPUTE_EXACT(lambda)
%   [ocv,gcv,fhat_r,f_r] = fcv.COMPUTE_EXACT(lambda)
%
% Input:
%   lambda - regularization parameter 
%
% Output:
%   ocv     - exact ordinary cross-validation score
%   gcv     - exact generalized cross-validation score
%   fhat_r  - corresponding Fourier coefficients
%   f_r     - corresponding function values

  [fhat_r,~] = lsqr(...
  @(x,transp_flag) A(self.plan,x,lambda,self.W,self.What,transp_flag),...
  [sqrt(self.W).*self.f;zeros(length(self.What),1)],1e-10);
  
  self.plan.fhat = fhat_r;
  nfsoft_trafo(self.plan);
  f_r = self.plan.f;
  
  h = zeros(self.M,1);
  for l = 1:self.M
    tmp = double( 1:self.M == l )';
    
    [fhat_tmp,~] = lsqr(...
    @(x,transp_flag) A(self.plan,x,lambda,self.W,self.What,transp_flag),...
    [sqrt(self.W).*tmp;zeros(length(self.What),1)]);
    
    self.plan.fhat = fhat_tmp;
    nfsoft_trafo(self.plan);
    tmp = self.plan.f;
    
    h(l) = tmp(l);
  end
  
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
