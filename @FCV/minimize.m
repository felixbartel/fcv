function lambda_min = minimize(self)
% fcv.MINIMIZE minimizes the cross-validation score
%
% Syntax:
%   lambda_min = fcv.MINIMIZE()
%
% Output:
%   lambda_min - lambda minimizing the cross-validation score

lambda0 = 1e-3;

options = optimoptions('fminunc',...
  'SpecifyObjectiveGradient',true,...
  'Algorithm','quasi-newton',...
  'Display','iter');
lambda_min = fminunc(@(t) objfcn(t,self),...
  log(lambda0),options);
lambda_min = exp(lambda_min);

end


function [val,grad] = objfcn(t,fcv)
  if nargout == 1
    s = fcv.compute(exp(t));
    val = s.gcv;
  else
    s = fcv.compute(exp(t),'derivative');
    val = s.gcv;
    grad = s.ddgcv;
    grad = real(grad/val*exp(t));
  end
  val = real(log(val));
end



