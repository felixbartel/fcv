function [f_r,fhat_r,ddf_r] = H(self,lambda)
% fcv.H computes the matrix-vector product with the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  fhat_r = alfft(self.z,self.M,self.I,self.W*self.f);
  fhat_r = fhat_r./(1+lambda*self.What);
  
  f_r = lfft(self.z,self.M,self.I,fhat_r.');
  if nargout > 2
    ddf_r = lfft(fhat_r.*self.What./(1+lambda*self.What).^2); % TODO: is this right?
  end
end
