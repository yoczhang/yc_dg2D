function getPolymeshOnLshape()
close all; clear; clc

Nelem = 1600;
NmaxIt = 1000;
meshName = ['Lshape_poly_',num2str(64)];

%% myLshapeDomain: [-1,1, -1,1]\[0,1, 0,1]
[node,elem,~,~,~]= PolyMesher(@myLshapeDomain,Nelem,NmaxIt);
% save(meshName,'node','elem')

%% fixed node
% y=-1
bool_idx = abs(node(:,2) - (-1))<1e-2;
node(bool_idx,2) = -1;

% y=1
bool_idx = abs(node(:,2) - (1))<1e-2;
node(bool_idx,2) = 1;

% x=-1
bool_idx = abs(node(:,1) - (-1))<1e-2;
node(bool_idx,1) = -1;

% x=1
bool_idx = abs(node(:,1) - (1))<1e-2;
node(bool_idx,1) = 1;

% x=0 && y>0
bool_idx1 = abs(node(:,1) - (0))<1e-2;
bool_idx2 = node(:,2) > 0;
bool_idx = bool_idx1 & bool_idx2;
node(bool_idx,1) = 0;

% y=0 && x>0
bool_idx1 = abs(node(:,2) - (0))<1e-2;
bool_idx2 = node(:,1) > 0;
bool_idx = bool_idx1 & bool_idx2;
node(bool_idx,2) = 0;

% additional condition
bool_idx1 = abs(node(:,2) - (1))<1e-6;
bool_idx2 = node(:,1) > 0;
bool_idx = bool_idx1 & bool_idx2;
node(bool_idx,1) = 0;

bool_idx1 = abs(node(:,1) - (1))<1e-6;
bool_idx2 = node(:,2) > 0;
bool_idx = bool_idx1 & bool_idx2;
node(bool_idx,2) = 0;

%%
patchPlotMesh(node,elem)
save(meshName,'node','elem')

end