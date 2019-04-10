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

  [fhat_r,~] = lsqr(...
    @(x,transp_flag) A(self.plan,x,lambda,self.W,self.What,transp_flag),...
    [sqrt(self.W).*self.f;zeros(length(self.What),1)],1e-10);
  
  nfsftmex('set_f_hat_linear',self.plan,fhat_r);
  nfsftmex('trafo',self.plan);
  f_r = nfsftmex('get_f',self.plan);
  
  h = sum((2*(0:self.N)'+1)./(1+lambda*self.What((1:self.N+1).^2)))/(4*pi)*self.W;
  ocv = real(sum(((f_r-self.f)./(1-h)).^2));
  gcv = real(sum(((f_r-self.f)./(1-mean(h))).^2));
end


function y = A(plan,x,lambda,W,W_hat,transp_flag)
  if strcmp(transp_flag,'notransp')
    nfsftmex('set_f_hat_linear',plan,x);
    nfsftmex('trafo',plan);
    y = nfsftmex('get_f',plan);
    
    y = sqrt(W).*y;
    
    y = [y; sqrt(lambda*W_hat).*x];
  else
    y = sqrt(W).*x(1:length(W));
    
    nfsftmex('set_f',plan,y);
    nfsftmex('adjoint',plan);
    y = nfsftmex('get_f_hat_linear',plan);
    
    y = y+sqrt(lambda*W_hat).*x(length(W)+1:end);
  end
end
