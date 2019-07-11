function [f_r,fhat_r,ddf_r] = H(self,lambda)
% fcv.H computes the matrix-vector product with the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  nfsftmex('set_f',self.plan,self.W.*self.f);
  nfsftmex('adjoint',self.plan);
  fhat_r = nfsftmex('get_f_hat_linear',self.plan);

  fhat_r = fhat_r./(1+lambda*self.What);
  nfsftmex('set_f_hat_linear',self.plan,fhat_r);
  nfsftmex('trafo',self.plan);
  f_r = nfsftmex('get_f',self.plan);
  if nargout > 2
    tmp = -fhat_r.*self.What./(1+lambda*self.What).^2;
    nfsftmex('set_f_hat_linear',self.plan,tmp);
    nfsftmex('trafo',self.plan);
    ddf_r = nfsftmex('get_f',self.plan);
  end
end
