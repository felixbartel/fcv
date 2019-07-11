function [h,ddh] = diagonals(self,lambda,~)
  h = sum(1./(1+lambda*self.What))*self.W;
  if nargout > 1
    ddh = self.W*sum(self.What./(1+lambda*self.What).^2);
  end
end
