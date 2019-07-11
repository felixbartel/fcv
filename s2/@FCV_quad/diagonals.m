function [h,ddh] = diagonals(self,lambda,~)
% fcv.DIAGONALS computes the diagonals of the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  h = sum((2*(0:self.N)'+1)./(1+lambda*self.What((1:self.N+1).^2)))/(4*pi)*self.W;
  if nargout > 1
    ddh = -sum((2*(0:self.N)'+1)./(1+lambda*self.What((1:self.N+1).^2).^2))/(4*pi)*self.W;
  end
end
