classdef FCV_appr < FCV
% FCV_APPR is a class for fast cross-validation without given quadrature rule
% on the rotation group
%
% Syntax:
%   fcv = FCV_APPR(nodes,f,W,N,s)
%   fcv = FCV_APPR(nodes,f,W,What)
%
% Input:
%   nodes - nodes in space domain
%   f     - function values
%   W     - weights in space domain
%   N     - bandwidth
%   s     - decay of Fourier coefficients
%   What  - weights in Fourier domain

properties
  nodes = []; % nodes in space domain
  f = [];     % function values
  W = [];     % weights in space domain
  What = [];  % weights in frequency domain
  plan = [];  % nfft plan
end

properties(Dependent = true)
  M;  % number of nodes
  N;  % bandwidth
end

methods
  function self = FCV_appr(angles,f,W,N,s)
    self.nodes = angles;
    self.f = f;
    if length(W) == 1
      self.W = W*ones(size(f));
    else
      self.W = W;
    end
    if length(N) == 1 % only bandwidth and decay are given
      self.What = ((0:N)+0.5).^(2*s);
      self.What = repelem(self.What,(1:2:(2*N+1)).^2)';
    else % all What are given
      self.What = N;
    end
    
    % parameters for compatibility to mtex
    self.plan = nfsoft(self.N,self.M,bitor(2^4,4),0,4,1000,2*ceil(1.5*self.N));
    
    self.plan.x = angles.';
  end

  function M = get.M(self)
    M = size(self.nodes,1);
  end
  
  function N = get.N(self)
    r = roots([4/3,4,11/3,1-length(self.What)]);
    N = r(real(r)-r == 0);
    N = round(N);
  end
end
 
end
