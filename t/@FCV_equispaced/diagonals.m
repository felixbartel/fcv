function [h,ddh] = diagonals(self,lambda,~)
% fcv.DIAGONALS computes the diagonals of the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  h = sum(1./(1+lambda*self.What))*self.W;
  if nargout > 1
    ddh = self.W*sum(self.What./(1+lambda*self.What).^2);
  end
end
