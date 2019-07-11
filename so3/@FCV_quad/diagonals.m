function [h,ddh] = diagonals(self,lambda,~)
% fcv.DIAGONALS computes the diagonals of the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  tmp = 2*(0:self.N)+1;
  tmp = repelem(tmp,(1:2:(2*self.N+1)).^2)';

  idx = nfsoft_f_hat_size(0:self.N);
  h = sum(tmp(idx)./(8*pi^2./tmp(idx)+lambda*self.What(idx)))*self.W;
%  if nargout > 1
%    ddh = self.W*sum(self.What./(1+lambda*self.What).^2);
%  end
end
