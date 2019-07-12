function [f_r,fhat_r,ddf_r] = H(self,lambda)
% fcv.H computes the matrix-vector product with the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  [fhat_r,~] = lsqr(...
    @(x,transp_flag) Afun(self.plan,x,lambda,self.W,self.What,transp_flag),...
    [sqrt(self.W).*self.f;zeros(self.M,1)],1e-10,200);
  
  tmp = fhat_r;
  tmp(1) = sqrt(2)*tmp(1);
  f_r = sqrt(self.N/2)*ndctIII(self.plan,tmp);
  
  if nargout > 2
% the derivative runs into numerical issues
    tmp = self.What.*fhat_r;
    [fhat_r,~] = lsqr(...
      @(x,transp_flag) Afun(self.plan,x,lambda,self.W,self.What,transp_flag),...
      [zeros(length(self.W),1);sqrt(self.What/lambda).*fhat_r],1e-10,200);
    
    tmp(1) = sqrt(2)*tmp(1);
    ddf_r = -sqrt(self.N/2)*ndctIII(self.plan,tmp);  
  end
end

function y = Afun(plan,x,lambda,W,What,transp_flag)
  M = length(W);
  N = length(What);
  if strcmp(transp_flag,'notransp')
    x(1) = sqrt(2)*x(1);
    y = sqrt(N/2*W).*ndctIII(plan,x);
    y = [y; sqrt(lambda*What).*x];
  else
    y = sqrt(N/2)*ndctII(plan,sqrt(W).*x(1:M));
    y(1) = sqrt(2)*y(1);
    y = y+sqrt(lambda*What).*x(M+1:end);
  end
end
