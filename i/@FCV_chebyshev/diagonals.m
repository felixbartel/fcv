function [h,ddh] = diagonals(self,lambda,~)
% fcv.DIAGONALS computes the diagonals of the hat matrix F*inv(F'*W*F+lambda*What)*F'*W
  b = 1./(pi*[1; 1/2*ones(self.N-1,1)]+lambda*self.What);

  btilde = zeros(2*self.M+1,1);
  btilde(1:2:2*self.M-1) = b;
  btilde(1) = btilde(1)*2*sqrt(2);

  h = sqrt(self.N/2)*dctI(btilde);
  h = h(2:2:2*self.M);
  h = h+(sum(b(2:end)));
  
  h = self.W/2*h;
 % if nargout > 1
 %   ddh = self.W*sum(self.What./(1+lambda*self.What).^2);
 % end
end
