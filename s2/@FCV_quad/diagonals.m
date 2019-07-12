function [h,ddh] = diagonals(self,lambda,~)
% fcv.DIAGONALS computes the diagonals of the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  idx = (1:self.N+1).^2;
  h = sum((2*(0:self.N)'+1)./(1+lambda*self.What(idx)))/(4*pi)*self.W;
  if nargout > 1
    ddh = -sum(self.What(idx).*(2*(0:self.N)'+1)./(1+lambda*self.What(idx)).^2)/(4*pi)*self.W;
  end
end
