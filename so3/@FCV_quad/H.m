function [f_r,fhat_r,ddf_r] = H(self,lambda)
% fcv.H computes the matrix-vector product with the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  self.plan.f = self.W.*self.f;
  nfsoft_adjoint(self.plan);
  fhat_r = self.plan.fhat;

  tmp = 2*(0:self.N)+1;
  tmp = repelem(tmp,(1:2:(2*self.N+1)).^2)';
  
  fhat_r = fhat_r./(8*pi^2./tmp+lambda*self.What);
  self.plan.fhat = fhat_r;
  nfsoft_trafo(self.plan);
  f_r = self.plan.f;
%  if nargout > 2
%    ddf_r = fftd(fhat_r.*self.What./(1+lambda*self.What).^2,self.d);
%  end
end
