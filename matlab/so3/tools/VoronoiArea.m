function W = VoronoiArea(nodes)
% VORONOIAREA computes the Voronoi weights for given nodes on SO(3)
%
% Syntax
%   W = VORONOIAREA(nodes)
%
% Input
%   nodes - M x d vecor
%
% Output
%   W - Voronoi weights

nodes = quaternion(reshape(double(nodes),[],4)');

M = length(nodes);
[v,c] = nodes.calcVoronoi;

W = zeros(M,1);
for idx = 1:M
  b1 = reshape(double(nodes(idx)),[],1);
  b1 = b1/norm(b1);
  
  b2 = reshape(double(v(c{idx}(1))),[],1)-b1;
  b2 = b2...
    -(b1.'*b2)*b1;
  b2 = b2/norm(b2);
  
  b3 = reshape(double(v(c{idx}(2))),[],1)-b1;
  b3 = b3...
    -(b1.'*b3)*b1...
    -(b2.'*b3)*b2;
  b3 = b3/norm(b3);
  
  b4 = reshape(double(v(c{idx}(3))),[],1)-b1;
  b4 = b4...
    -(b1.'*b4)*b1...
    -(b2.'*b4)*b2...
    -(b3.'*b4)*b3;
  b4 = b4/norm(b4);
  
  V = reshape(double(v(c{idx})),[],4)';
  V = (([b1 b2 b3 b4])\V)'; 
  V = diag(1-2*( V(:,1) < 0 ))*V;
  V = V(:,2:end);
  [~,W(idx)] = convhulln(V);
end

W = 8*W;
  
end


% %% S2
% 
% nodes = vector3d.rand(1000);
% 
% [v,c] = nodes.calcVoronoi;
% 
% s = zeros(length(nodes),1);
% for idx = 1:length(nodes)
%   b1 = nodes(idx).xyz;
%   b1 = b1/norm(b1);
%   
%   b2 = v(c{idx}(1)).xyz-b1;
%   b2 = b2...
%     -(b1*b2.')*b1;
%   b2 = b2/norm(b2);
%   
%   b3 = v(c{idx}(2)).xyz-b1;
%   b3 = b3...
%     -(b1*b3.')*b1...
%     -(b2*b3.')*b2;
%   b3 = b3/norm(b3);
%   
%   
%   V = v(c{idx});
%   V = (([b1' b2' b3'])\V.xyz')';
%   V = V(:,2:end);
%   [~,s(idx)] = convhulln(V);
% end
% abs(sum(s)-4*pi)/(4*pi)
% s = nodes.calcVoronoiArea();
% abs(sum(s)-4*pi)/(4*pi)




