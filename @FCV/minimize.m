function lambda_min = minimize(self)
% fcv.MINIMIZE minimizes the cross-validation score
%
% Syntax:
%   lambda_min = fcv.MINIMIZE()
%
% Output:
%   lambda_min - lambda minimizing the cross-validation score

lambda0 = 0.1;

if ismethod(self,'compute_with_grad')
  MaxFunEvals = 7; % because seven is a lucky number
  options = optimoptions('fminunc',...
    'SpecifyObjectiveGradient',true,...
    'Algorithm','trust-region',...
    'Display','iter',...
    'MaxFunEvals',MaxFunEvals);
else
  MaxFunEvals = 50;
  options = optimoptions('fminunc',...
    'SpecifyObjectiveGradient',false,...
    'Algorithm','quasi-newton',...
    'Display','iter',...
    'MaxFunEvals',MaxFunEvals);
end
lambda_min = fminunc(@(t) objfcn(t,self),...
  log(lambda0), options);
lambda_min = exp(lambda_min);

end


function [val,grad] = objfcn(t,fcv)
  if nargout == 1
    s = fcv.compute(exp(t));
    val = s.gcv;
  else
    s = fcv.compute_with_grad(exp(t));
    val = s.gcv;
    grad = s.gcv_grad;
    grad = real(grad/val*exp(t));
  end
  val = real(log(val));
end



