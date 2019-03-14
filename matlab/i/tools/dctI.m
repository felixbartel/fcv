function ahat = dctI(a)
  s = size(a);
  a = a(:);
  N = length(a);

  y = a([1:N N-1:-1:2]);
  y([1 N]) = sqrt(2)*y([1 N]);

  yhat = fft(y);

  ahat = 1/sqrt(2*N-2)*real(yhat(1:N));
  ahat([1 end]) = ahat([1 end])/sqrt(2);

  reshape(ahat,s);
end
