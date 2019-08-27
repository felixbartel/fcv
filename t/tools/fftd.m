function f = fftd(fhat,d)
% d-dimensional fft with linear in- and output
  N = (length(fhat))^(1/d);
  if d == 1
    fhat = fhat([N/2+1:end 1:N/2]);
    f = fft(fhat);
  else
    f = reshape(fhat,N*ones(1,d));
    
    for j = 1:d
      f = f([N/2+1:end 1:N/2],:);
      f = reshape(f,N*ones(1,d));
      f = shiftdim(f,1);
    end
    
    for j = 1:d
      f = fft(f,[],j);
    end
    f = reshape(f,[],1);
  end
end
