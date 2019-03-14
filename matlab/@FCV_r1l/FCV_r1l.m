classdef FCV_r1l
% FCV_r1L is a class for fast cross-validation for rank-1 lattices
% on the torus
%
% Syntax:
%   fcv = FCV_R1L(z,M,f,I,What)
%
% Input:
%   z    - generating vector
%   M    - number of nodes
%   f    - function values
%   I    - Index set
%   What - weights in Fourier space

properties
  z = [];     % generating vector
  M = [];     % number of nodes
  f = [];     % function values
  W = [];     % weights in space domain
  I = [];     % index set in frequency domain
  What = [];  % weights in frequency domain
end

properties(Dependent = true)
  d;  % dimension
end

methods
  function self = FCV_r1l(z,M,f,I,What)
    self.z = z;
    self.M = M;
    self.f = f;
    self.W = 1/M;
    self.I = I;
    self.What = What;
  end
  
  function d = get.d(self)
    d = length(self.z);
  end
end
 
end
