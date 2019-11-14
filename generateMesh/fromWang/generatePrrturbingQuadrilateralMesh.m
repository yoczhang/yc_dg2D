function [node, elem] = generatePrrturbingQuadrilateralMesh(node, elem, area, theta)

N = size(node,1);
bdNodeIdx = findboundary(elem);

% Note that the interface location is y = 1.
[idx, unused] = find(node(:,2)==1);
bdNodeIdx = union(bdNodeIdx, idx);

patchArea = sparse(elem,ones(size(elem,1),4),area*[1,1,1,1],N,1);
h = sqrt(patchArea/6);
dp = theta*[h h].*(2*rand(N,2)-1);
oldNode = node(bdNodeIdx,:);
node = node + dp;
node(bdNodeIdx,:) = oldNode;

end