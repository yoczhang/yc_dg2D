function [node,elem] = triMeshOnPolygonDomain(X1, X2, h)
%
%   This function allows to generate a triangulation for arbitrary polygonal
%	domains.
%
%   input: 
%       X1,X2, the (x^1,x^2) coordinates of the boundary points 
%                 that describe the domain.
%       h, the maximum diameter of an element.
%
%   output:
%       node, the (x,y)-coordinate of mesh nodes.
%       elem, 
%
%   YcZhang 7/10/2017
%
%   Last modified 7/10/2017
%

assert(length(X1) >= 3, 'At least 3 points are required for a 2D domain')
assert(isequal(size(X1), size(X2)), 'X1 and X2 must be of same size')
gd = [2; length(X1(:)); X1(:); X2(:)]; % geometry description
sf = 'polygon';                        % set formula
ns = double('polygon')';               % name space
[p, e, t] = initmesh(decsg(gd,sf,ns), 'Hmax', h);
node = p';
elem = t(1:3,:)';

end % function