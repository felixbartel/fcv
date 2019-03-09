function y = fft2(y,flag)
% two-dimensional fft
  N = sqrt(length(y));
  if strcmp(flag,'notransp')
    y = reshape(fft(fft(reshape(y,N,N),[],2)),[],1);
  elseif strcmp(flag,'transp')
    y = N^2*reshape(ifft(ifft(reshape(y,N,N),[],2)),[],1);
  end
end