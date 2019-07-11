function h = diagonals(self,lambda,flag)
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
    tmp = 2*(0:self.N)+1;
    tmp = repelem(tmp,(1:2:(2*self.N+1)).^2)';

    idx = nfsoft_f_hat_size(0:self.N);
    h = sum(tmp(idx)./(8*pi^2./tmp(idx)+lambda*self.What(idx)))*self.W;
  end
end
