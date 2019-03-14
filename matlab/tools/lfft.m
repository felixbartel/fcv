function f = lfft(z,M,I,fhat)
% LFFT rank-1 lattice based Fourier transformation
%
% Syntax:
%   f = LFFT(z,M,I,fhat)
%
% Input:
%   z     - generating vector
%   M     - number of nodes
%   I     - Index set
%   fhat  - Fourier coefficients
%
% Output:
%   f - function values

  k = mod(I*z.',M)+1;
  k = floor(k);
  fhat1 = accumarray(k,fhat,[M,1],@sum);
  f = M*ifft(fhat1);
end
