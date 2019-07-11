function [f_r,fhat_r,ddf_r] = H(self,lambda)
% fcv.H computes the matrix-vector product with the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  [fhat_r,~] = lsqr(...
    @(x,transp_flag) Afun(self.plan,x,lambda,self.W,self.What,transp_flag),...
    [sqrt(self.W).*self.f;zeros(length(self.What),1)],1e-10);
    
  nfft_set_f_hat(self.plan,fhat_r);
  nfft_trafo(self.plan);
  f_r = nfft_get_f(self.plan);

  if nargout > 2
    [tmp,~] = lsqr(...
      @(x,transp_flag) Afun(self.plan,x,lambda,self.W,self.What,transp_flag),...
      [zeros(length(self.W),1);-self.What.*fhat_r],1e-10);
    nfft_set_f_hat(self.plan,tmp);
    nfft_trafo(self.plan);
    ddf_r = nfft_get_f(self.plan);
  end
end

function y = Afun(plan,x,lambda,W,What,transp_flag)
  if strcmp(transp_flag,'notransp')
    nfft_set_f_hat(plan,x);
    nfft_trafo(plan);
    y = nfft_get_f(plan);
    
    y = sqrt(W).*y;
    
    y = [y; sqrt(lambda*What).*x];
  else
    y = sqrt(W).*x(1:length(W));
    
    nfft_set_f(plan,y);
    nfft_adjoint(plan);
    y = nfft_get_f_hat(plan);
    
    y = y+sqrt(lambda*What).*x(length(W)+1:end);
  end
end
