function [f_r,fhat_r,ddf_r] = H(self,lambda)
% fcv.H computes the matrix-vector product with the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  fhat_r = dctII(self.W*self.f);
  fhat_r(1) = sqrt(2)*fhat_r(1);
  fhat_r = sqrt(self.N/2)*fhat_r;
  
  fhat_r = fhat_r./(pi*[1; 1/2*ones(self.N-1,1)]+lambda*self.What);
  
  fhat_r(1) = sqrt(2)*fhat_r(1);
  f_r = sqrt(self.N/2)*dctIII(fhat_r);
 % if nargout > 2
 %   ddf_r = fftd(fhat_r.*self.What./(1+lambda*self.What).^2,self.d);
 % end
end
