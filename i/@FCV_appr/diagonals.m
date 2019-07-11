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
  % compute diagonal emelents
    b = 1./(pi*[1; 1/2*ones(self.N-1,1)]+lambda*self.What);
    b(1) = 2*b(1);

    bhat_r = zeros(2*self.M,1);
    bhat_r(1:2:2*self.M-1) = b;

    nfct_set_f_hat(self.plan2,double(bhat_r));
    nfct_trafo(self.plan2);
    h = nfct_get_f(self.plan2);

    h = self.W.*(h+sum(b(2:end)))/2;
  end
end
