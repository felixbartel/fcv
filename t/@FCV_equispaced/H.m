function [f_r,fhat_r,ddf_r] = H(self,lambda)
% fcv.H computes the matrix-vector product with the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  fhat_r = fftd_adj(self.W*self.f,self.d);
  fhat_r = fhat_r./(1+lambda*self.What);
  f_r = fftd(fhat_r,self.d);
  if nargout > 2
    ddf_r = -fftd(fhat_r.*self.What./(1+lambda*self.What),self.d);
  end
end
