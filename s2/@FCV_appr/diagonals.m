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
    idx = (1:self.N+1).^2;
    h = sum((2*(0:self.N)'+1)./(1+lambda*self.What(idx)))/(4*pi)*self.W;
    if nargout > 1
      ddh = -sum(self.What(idx).*(2*(0:self.N)'+1)./(1+lambda*self.What(idx)).^2)/(4*pi)*self.W;
    end
  end
end
