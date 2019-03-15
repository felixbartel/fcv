function [ocv,gcv,fhat_r,f_r] = compute_exact(self,lambda)
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
    [sqrt(self.M/2*self.W).*self.f;zeros(self.M,1)],1e-8,1000);
  
  f_r = ndctIII(self.plan,fhat_r); 
   
  h = zeros(self.M,1);
  for l = 1:self.M
    tmp = double( 1:self.M == l )';
    
    [fhat_tmp,~] = lsqr(...
    @(x,transp_flag) A(self.plan,x,lambda,self.W,self.What,transp_flag),...
    [sqrt(self.M/2*self.W).*tmp;zeros(self.M,1)],1e-8,1000);
    
    tmp = ndctIII(self.plan,fhat_tmp);
    
    h(l) = tmp(l);
  end

  ocv = real(sum(((f_r-self.f)./(1-h)).^2));
  gcv = real(sum(((f_r-self.f)./(1-mean(h))).^2));
end


function y = A(plan,x,lambda,W,What,transp_flag)
  M = length(W);
  if strcmp(transp_flag,'notransp')
    y = sqrt(M/2*W).*ndctIII(plan,x);
    y = [y; sqrt(lambda*What).*x];
  else
    y = ndctII(plan,sqrt(M/2*W).*x(1:M));
    y = y+sqrt(lambda*What).*x(M+1:end);
  end
end
