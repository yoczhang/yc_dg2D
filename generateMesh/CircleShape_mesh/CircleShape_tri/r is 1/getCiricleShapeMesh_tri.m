function getCiricleShapeMesh_tri()
%
close all; clear; clc

hh = 64;
initH = 1/hh;
meshName = ['Circleshape_tri_',num2str(hh)];

[node,elem] = circlemesh(0,0,1,initH);
%%
patchPlotMesh(node,elem)
save(meshName,'node','elem')
end