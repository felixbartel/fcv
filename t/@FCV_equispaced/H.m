function [f_r,fhat_r,ddf_r] = H(self,lambda)
  fhat_r = self.W*fftd_adj(self.f,self.d);
  fhat_r = fhat_r./(1+lambda*self.What);
  f_r = fftd(fhat_r,self.d);
  if nargout > 2
    ddf_r = fftd(fhat_r.*self.What./(1+lambda*self.What).^2,self.d);
  end
end
