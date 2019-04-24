function f_hat = ndctII(plan,f)
  M = length(f);
  nfct_set_f(plan,double(f));
  nfct_adjoint(plan);
  f_hat_adjoint = nfct_get_f_hat(plan);

  f_hat = sqrt(2/M)*f_hat_adjoint;
  f_hat(1) = f_hat(1)/sqrt(2);
end
