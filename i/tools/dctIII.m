function ahat = dctIII(a)
  s = size(a);
  a = a(:);
  N = length(a);

  y = [sqrt(2)*a(1); ...
    exp(-2i*pi*(1:N-1)'/(4*N)).*a(2:N); ...
    0; ...
    exp(-2i*pi*((3*N+1):(4*N-1))'/(4*N)).*a(N:-1:2)];
  yhat = fft(y);
  ahat = 1/sqrt(2*N)*real(yhat(1:N));

  ahat = reshape(ahat,s);
end
