function getCircleShapeMesh_poly()
close all; clear; clc

Nelem = 1600;
NmaxIt = 900;
meshName = ['Circleshape_poly_',num2str(64)];

%% myCircleDomain: x^2+y^2 <= 1
[node,elem,~,~,~]= PolyMesher(@myCircleShape,Nelem,NmaxIt);

%%
patchPlotMesh(node,elem)
save(meshName,'node','elem')
end