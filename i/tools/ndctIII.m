function f = ndctIII(plan,f_hat)
  M = length(f_hat);
  f_hat = [f_hat(1)/sqrt(2); f_hat(2:end)];
  nfct_set_f_hat(plan,double(f_hat));
  nfct_trafo(plan);
  f = nfct_get_f(plan);

  f = sqrt(2/M)*f;
end
