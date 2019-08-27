function fhat = fftd_adj(f,d)
% d-dimensional adjoint fft with linear in- and output
  N = (length(f))^(1/d);
  if d == 1
    fhat = N*ifft(f);
    fhat = fhat([N/2+1:end 1:N/2]);
  else
    fhat = reshape(f,N*ones(1,d));
    for j = 1:d
      fhat = N*ifft(fhat,[],j);
    end
    for j = 1:d
      fhat = fhat([N/2+1:end 1:N/2],:);
      fhat = reshape(fhat,N*ones(1,d));
      fhat = shiftdim(fhat,1);
    end
    fhat = reshape(fhat,[],1);
  end
end
