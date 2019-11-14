%% Illustrations of all used meshes
% illustrationOfMeshes;
clearvars;
close all
clc
%----------- yc test -----------------
%------ mesh1 
%-------------------- Quad mesh ---------------------
n = 1;
h_x = 1/2^(n+1); h_y = h_x;
h_partition=[h_x,h_y];
left=0;right=1;bottom=0;top=1;
[M,T]=generate_quad_P_T(left,right,bottom,top,h_partition,1);
node1 = M'; elem1 = T';

[node,elem] = squarequadmesh([0,1,0,1],1/2^2); 
[nodeNew,elemNew] = non_convex_octagona_mesh(node,elem);




%----------------------------------------------------
nodeNew;
elemNew;

