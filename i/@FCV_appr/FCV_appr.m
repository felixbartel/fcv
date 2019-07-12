classdef FCV_appr < FCV
% FCV_APPR is a class for fast cross-validation for nonequispaced
% nodes on the unit interval
%
% Syntax:
%   fcv = FCV_APPR(nodes,f,W,N,s)
%   fcv = FCV_APPR(nodes,f,W,What)
%
% Input:
%   nodes - nodes in space domain
%   f     - function values
%   W     - weights in space domain (if empty, adopted VoronoiArea is used)
%   N     - bandwidth
%   s     - decay of Fourier coefficients
%   What  - weights in Fourier domain

properties
  nodes = []; % nodes in space domain
  f = [];     % function values
  W = [];     % weights in space domain
  What = [];  % weights in frequency domain
  plan = [];  % nfct plan 1
  plan2 = []; % nfct plan 2
  x0 = [];    % initial guess while solving with H
end

properties(Dependent = true)
  M;  % number of nodes
  N;  % bandwidth
end

methods
  function self = FCV_appr(nodes,f,W,N,s)
    self.nodes = nodes;
    self.f = f;
    if isempty(W)
      [nodes_tilde,idx] = sort(nodes);
      nodes_tilde = acos(nodes_tilde);

      self.W = pi-(nodes_tilde(1)+nodes_tilde(2))/2;                    % first node
      self.W = [self.W; (nodes_tilde(1:end-2)-nodes_tilde(3:end))/2];   % indermediate nodes
      self.W = [self.W; (nodes_tilde(end-1)+nodes_tilde(end))/2];       % last nodes
      self.W(idx) = self.W;
    else
      self.W = W;
    end
    
    
    if length(N) == 1 % only bandwidth and decay are given
      self.What  = ((1:N).').^s;
    else % all What are given
      self.What = What;
    end
    
    self.plan = nfct_init_1d(self.N,self.M);
    nfct_set_x(self.plan,acos(self.nodes.')/(2*pi));

    self.plan2 = nfct_init_1d(2*self.N,self.M);
    nfct_set_x(self.plan2,acos(self.nodes.')/(2*pi));
  end
    
  function delete(self)
    nfct_finalize(self.plan);
    nfct_finalize(self.plan2);
  end
  
  function M = get.M(self)
    M = length(self.f);
  end
  
  function N = get.N(self)
    N = length(self.What);
  end
end
 
end
