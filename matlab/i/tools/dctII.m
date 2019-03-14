function ahat = dctII(a)
  s = size(a);
  a = a(:);
  N = length(a);

  y = a([1:N N:-1:1]);

  yhat = fft(y);

  ahat = 1/sqrt(2*N)*real(exp(-2i*pi*(0:N-1)'/(4*N)).*yhat(1:N));
  ahat(1) = ahat(1)/sqrt(2);

  reshape(ahat,s);
end
