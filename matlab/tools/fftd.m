function fhat = fftd(f,d)
% d-dimensional fft with linear in- and output
  if d == 1
    fhat = fft(f);
  else
    N = (length(f))^(1/d);
    fhat = reshape(f,N*ones(1,d));
    for j = 1:d
      fhat = fft(fhat,[],j);
    end
    fhat = reshape(fhat,[],1);
  end
end