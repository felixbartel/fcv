function lambda_min = minimize(fcv)
% fcv.MINIMIZE minimizes the cross-validation score
%
% Syntax:
%   lambda_min = fcv.MINIMIZE()
%
% Output:
%   lambda_min - lambda minimizing the cross-validation score

lambda0 = 1;
MaxFunEvals = 7; % because seven is a lucky number

options = optimoptions('fminunc',...
  'SpecifyObjectiveGradient',true,...
  'Algorithm','trust-region',...
  'Display','iter',...
  'MaxFunEvals',MaxFunEvals);

lambda_min = fminunc(@(t) objfcn(t,fcv),...
  log(lambda0), options);
lambda_min = exp(lambda_min);

end


function [val,grad] = objfcn(t,fcv)
  if nargout == 1
    val = fcv.compute(exp(t));
  else
    [val,grad] = fcv.compute_with_grad(exp(t));
    grad = grad/val*exp(t);
  end
  val = log(val);
end


