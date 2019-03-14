classdef FCV_chebyshev < FCV
% FCV_CHEBYSHEV is a class for fast cross-validation for Chebyshev nodes
% on the unit interval
%
% Syntax:
%   fcv = FCV_CHEBYSHEV(f,s)
%   fcv = FCV_CHEBYSHEV(f,What)
%
% Input:
%   f     - function values
%   s     - decay of Fourier coefficients
%   What  - weights in Fourier domain

properties
  f = [];     % function values
  W = [];     % weights in space domain
  What = [];  % weights in frequency domain
end

properties(Dependent = true)
  M;  % number of nodes
  N;  % bandwidth
end

methods
  function self = FCV_chebyshev(f,What)
    self.f = f;
    self.W = 2/length(f);
    if length(What) == 1 % only decay is given
      self.What  = ((1:self.M).').^What;
    else % all What are given
      self.What = What;
    end
  end
  
  function M = get.M(self)
    M = length(self.f);
  end
  
  function N = get.N(self)
    N = length(self.What);
  end
end
 
end
