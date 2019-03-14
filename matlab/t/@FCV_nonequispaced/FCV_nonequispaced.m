classdef FCV_nonequispaced < FCV
% FCV_NONEQUISPACED is a class for fast cross-validation for equispaced nodes
% on the torus
%
% Syntax:
%   fcv = FCV_NONEQUISPACED(nodes,f,W,N,s)
%   fcv = FCV_NONEQUISPACED(nodes,f,W,What)
%
% Input:
%   nodes - nodes in space domain
%   f     - function values
%   W     - weights in space domain (if empty, VoronoiArea is used)
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
  d;  % dimension
  M;  % number of nodes
  N;  % bandwidth
end

methods
  function self = FCV_nonequispaced(nodes,f,W,N,s)
    self.nodes = nodes;
    self.f = f;
    if isempty(W)
      self.W = VoronoiArea(nodes);
    else
      self.W = W;
    end
    if length(N) == 1 % only bandwidth and decay are given
      self.What = zeros(N*ones(1,self.d));
      t = [-N/2:N/2-1];
      if self.d == 1
        self.What = t.^2;
      else
        for j = 1:self.d
          t = reshape(t,(N-1)*((1:self.d) == j)+ones(1,self.d));
          self.What = self.What+t.^2;
        end
      end
      self.What = 1+self.What.^(s/2);
      self.What = self.What(:);
    else % all What are given
      self.What = N;
    end
    
    self.plan = nfft_init(self.N*ones(1,self.d),self.M);
    nfft_set_x(self.plan,self.nodes');
    nfft_precompute_psi(self.plan);  
  end
  
  function delete(self)
    nfft_finalize(self.plan);
  end
  
  function d = get.d(self)
    d = size(self.nodes,2);
  end
  
  function M = get.M(self)
    M = size(self.nodes,1);
  end
  
  function N = get.N(self)
    N = length(self.What)^(1/self.d);
  end
end
 
end
