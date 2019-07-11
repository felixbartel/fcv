function [h,ddh] = diagonals(self,lambda,~)
% fcv.DIAGONALS computes the diagonals of the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  h = sum((2*(0:self.N)'+1)./(1+lambda*self.What((1:self.N+1).^2)))/(4*pi)*self.W;
%  if nargout > 1
%    ddh = self.W*sum(self.What./(1+lambda*self.What).^2);
%  end
end
