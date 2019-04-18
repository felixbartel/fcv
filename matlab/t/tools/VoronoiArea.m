function W = VoronoiArea(nodes)
% VORONOIAREA computes the Voronoi weights for given nodes
%
% Syntax
%   W = VORONOIAREA(nodes)
%
% Input
%   nodes - M x d vecor (only d=1,2,3 supported so far)
%
% Output
%   W - Voronoi weights

d = size(nodes,2);
M = size(nodes,1);

if d == 1 % one-dimensional case
  [nodes,I] = sort(nodes);
  W = zeros(M,1);
  W(1) = nodes(2)-nodes(end)+1;
  W(2:end-1) = nodes(3:end)-nodes(1:end-2);
  W(end) = nodes(1)+1-nodes(end-1);
  W(I) = W/2;
  
else % d > 1
  % firstly make nodes peridic
  shift = -ones(1,d);
  nodes_p = nodes;
  while shift(1) < 2
    if norm(shift-zeros(1,d)) > 0
      nodes_p = [nodes_p; nodes+shift];
    end
    shift(end) = shift(end)+1;
    for idx = length(shift):-1:2
      if shift(idx) > 1
        shift(idx) = -1;
        shift(idx-1) = shift(idx-1)+1;
      end
    end
  end

  [V,C] = voronoin(nodes_p);
   
  W = zeros(M,1);
  for idx = 1:M % todo: does not work for d > 3
    [~,W(idx)] = convhull(V(C{idx},:));
  end
end
end
