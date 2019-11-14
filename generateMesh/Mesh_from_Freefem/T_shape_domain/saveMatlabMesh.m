%
%
%

[node, edge, elem]=importfilemesh('TshapeDomain1.msh');
node = node';
elem = elem(1:3,:)';

save TshapeDomain1 node elem