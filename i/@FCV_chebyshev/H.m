function [f_r,fhat_r,ddf_r] = H(self,lambda)
% fcv.H computes the matrix-vector product with the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  fhat_r = sqrt(self.N/2)*dctII(self.W*self.f);
  fhat_r(1) = sqrt(2)*fhat_r(1);
  
  fhat_r = fhat_r./(pi*[1; 1/2*ones(self.N-1,1)]+lambda*self.What);
  
  tmp = fhat_r;
  tmp(1) = sqrt(2)*tmp(1);
  f_r = sqrt(self.N/2)*dctIII(tmp);
  
  if nargout > 2
    tmp = self.What.*fhat_r./(pi*[1; 1/2*ones(self.N-1,1)]+lambda*self.What);

    tmp(1) = sqrt(2)*tmp(1);
    ddf_r = -sqrt(self.N/2)*dctIII(tmp);  
  end
end
