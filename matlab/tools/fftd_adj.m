function f = fftd_adj(fhat,d)
% d-dimensional adjoint fft with linear in- and output
  N = (length(fhat))^(1/d);
  if d == 1
    f = N*ifft(fhat);
  else
    f = reshape(fhat,N*ones(1,d));
    for j = 1:d
      f = N*ifft(f,[],j);
    end
    f = reshape(f,[],1);
  end
end