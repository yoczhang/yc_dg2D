function getCircleShapeMesh_poly()
close all; clear; clc

% in r = 0.5
% about (Nelem=50--mesh4),
% (Nelem=100--mesh8),
% (Nelem=200--mesh16), 
% (Nelem=400--mesh32), 
% (Nelem=800--mesh64), 
Nelem = 800;
NmaxIt = 800;
meshName = ['Circleshape_poly_',num2str(64)];

%% myCircleDomain: x^2+y^2 <= 1
[node,elem,~,~,~]= PolyMesher(@myCircleShape,Nelem,NmaxIt);

%%
patchPlotMesh(node,elem)
save(meshName,'node','elem')
end