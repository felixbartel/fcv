classdef FCV_equispaced < FCV
% FCV_EQUISPACED is a class for fast cross-validation for equispaced nodes
% on the torus
%
% Syntax:
%   fcv = FCV_EQUISPACED(d,f,s)
%   fcv = FCV_EQUISPACED(d,f,What)
%
% Input:
%   d     - dimension
%   f     - function values
%   s     - decay of Fourier coefficients
%   What  - weights in Fourier domain

properties
  d = [];     % dimension
  f = [];     % function values
  W = [];     % weights in space domain
  What = [];  % weights in frequency domain
end

properties(Dependent = true)
  M;  % number of nodes
  N;  % bandwidth
end

methods
  function self = FCV_equispaced(d,f,What)
    self.d = d;
    self.f = f;
    self.W = 1/length(f);
    if length(What) == 1 % only decay is given
      s = What;
      if d == 1
        self.What = (1+abs(-self.N/2:self.N/2-1).^s).';
      else
        self.What = zeros(self.N*ones(1,d));
        t = abs(-self.N/2:self.N/2-1);
        %t = [0:self.N/2-1 self.N/2:-1:1];
        for j = 1:d
          t = reshape(t,(self.N-1)*((1:d) == j)+ones(1,d));
          self.What = self.What+t;
        end
        self.What = 1+self.What.^s;
        self.What = self.What(:);
      end
    else % all What are given
      self.What = What;
    end
    
  end
  
  function M = get.M(self)
    M = length(self.f);
  end
  
  function N = get.N(self)
    N = length(self.f)^(1/self.d);
  end
end
 
end
