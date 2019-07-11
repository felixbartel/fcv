function [h,ddh] = diagonals(self,lambda,flag)
% fcv.DIAGONALS computes the diagonals of the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  if nargin < 3; flag = ""; end
  if strcmp(flag,"exact")
    f_copy = self.f;
    h = zeros(self.M,1);
    for l = 1:self.M
      self.f = double( 1:self.M == l )';
      tmp = self.H(lambda);
      h(l) = tmp(l);
    end
    self.f = f_copy;
  else
    h = sum(1./(1+lambda*self.What))*self.W;
    if nargout > 1
      ddh = -self.W*sum(self.What./(1+lambda*self.What).^2);
    end
  end
end
