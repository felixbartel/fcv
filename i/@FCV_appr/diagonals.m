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
    b = 1./(pi*[1; 1/2*ones(self.N-1,1)]+lambda*self.What);

    btilde = zeros(2*self.M+1,1);
    btilde(1:2:2*self.M-1) = b;
    btilde(1) = btilde(1)*2*sqrt(2);

    h = sqrt(self.N/2)*dctI(btilde);
    h = h(2:2:2*self.M);
    h(1) = h(1)*sqrt(2);
    h = h+(sum(b(2:end)));

    h = self.W/2.*h;

    if nargout > 1
      ddb = -self.What./(pi*[1; 1/2*ones(self.N-1,1)]+lambda*self.What).^2;
      ddbtilde = zeros(2*self.M+1,1);
      ddbtilde(1:2:2*self.M-1) = ddb;
      ddbtilde(1) = ddbtilde(1)*2*sqrt(2);

      ddh = sqrt(self.N/2)*dctI(ddbtilde);
      ddh = ddh(2:2:2*self.M);
      ddh(1) = ddh(1)*sqrt(2);
      ddh = ddh+(sum(ddb(2:end)));

      ddh = self.W/2.*ddh;
    end
end
