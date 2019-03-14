function fhat = alfft(z,M,I,f)
% ALFFT adjoined rank-1 lattice based Fourier transformation
%
% Syntax:
%   fhat = ALFFT(z,M,I,f)
%
% Input:
%   z     - generating vector
%   M     - number of nodes
%   I     - Index set
%   f     - function values
%
% Output:
%   fhat - Fourier coefficients
  k = mod(I*z.',M)+1;
  k = floor(k);
  ghat = fft(f);
  fhat = ghat(k);
end
