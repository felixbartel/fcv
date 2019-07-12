function [f_r,fhat_r,ddf_r] = H(self,lambda)
% fcv.H computes the matrix-vector product with the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  [fhat_r,~] = lsqr(...
    @(x,transp_flag) Afun(self.plan,x,lambda,self.W,self.What,transp_flag),...
    [sqrt(self.W).*self.f;zeros(length(self.What),1)],1e-10,[],[],[],self.x0);
  self.x0 = fhat_r;
  
  nfsftmex('set_f_hat_linear',self.plan,fhat_r);
  nfsftmex('trafo',self.plan);
  f_r = nfsftmex('get_f',self.plan);
  
  if nargout > 2
    [tmp,~] = lsqr(...
      @(x,transp_flag) Afun(self.plan,x,lambda,self.W,self.What,transp_flag),...
      [zeros(length(self.W),1);sqrt(self.What/lambda).*fhat_r],1e-10);
    nfsftmex('set_f_hat_linear',self.plan,tmp);
    nfsftmex('trafo',self.plan);
    ddf_r = -nfsftmex('get_f',self.plan);
  end
end

function y = Afun(plan,x,lambda,W,What,transp_flag)
  if strcmp(transp_flag,'notransp')
    nfsftmex('set_f_hat_linear',plan,x);
    nfsftmex('trafo',plan);
    y = nfsftmex('get_f',plan);
    
    y = sqrt(W).*y;
    
    y = [y; sqrt(lambda*What).*x];
  else
    y = sqrt(W).*x(1:length(W));
    
    nfsftmex('set_f',plan,y);
    nfsftmex('adjoint',plan);
    y = nfsftmex('get_f_hat_linear',plan);
    
    y = y+sqrt(lambda*What).*x(length(W)+1:end);
  end
end
