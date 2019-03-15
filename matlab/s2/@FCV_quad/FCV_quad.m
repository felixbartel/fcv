classdef FCV_quad < FCV
% FCV_QUAD is a class for fast cross-validation for given quadrature rule
% on the two-dimensional sphere
%
% Syntax:
%   fcv = FCV_QUAD(nodes,f,W,N,s)
%   fcv = FCV_QUAD(nodes,f,W,What)
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
  function self = FCV_quad(nodes,f,W,N,s)
    self.nodes = nodes;
    self.f = f;
    self.W = W;
    if length(N) == 1 % only bandwidth and decay are given
      self.What = (2*(0:N)+1).^(2*s);
      self.What = repelem(self.What,1:2:(2*N+1))';
    else % all What are given
      self.What = N;
    end
    
    nfsftmex('precompute',self.N,1000,1,0);
    self.plan = nfsftmex('init_advanced',self.N,self.M,1);
    nfsftmex('set_x',self.plan,[self.nodes(:,1)';self.nodes(:,2)']);
    nfsftmex('precompute_x',self.plan);
  end
  
  function delete(self)
    nfsft_finalize(self.plan);
  end

  function M = get.M(self)
    M = size(self.nodes,1);
  end
  
  function N = get.N(self)
    N = sqrt(length(self.What))-1;
  end
end
 
end
