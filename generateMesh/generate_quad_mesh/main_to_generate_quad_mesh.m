function main_to_generate_quad_mesh
%
%   Here we make the rules for naming the mesh, 
%   such as: 'S_quadmesh_4times2_[0_1]_[1-2_1]',
%       S_quadmesh, is the mesh name,
%       4time2, is the x-coordinate and y-coordinate elem arrange,
%       [0_1], is the x-coordinate domain [0,1]
%       [1-2_1], is the y-coordinate domain [1/2,1], and 1-2 is the transformation of 1/2.
%       
%
%   YcZhang 3/9/2017
%
%   Last modified 4/9/2017
%

clear
close all

% generate the Stokes-Darcy mesh, 
StokesDarcyDomain_1


end % function


%% sub function 1
function StokesDarcyDomain_1
% Stokes domain: [0,1]x[1/2,1]. Darcy domain: [0,1]x[0,1/2].
h_x = 1/64; h_y = h_x;
h=[h_x,h_y];

% Stokes domain
left_S=0;right_S=1;bottom_S=1/2;top_S=1;
[M_S,T_S]=generate_quad_P_T(left_S,right_S,bottom_S,top_S,h,1);
elements = T_S';
vertices = M_S';
save S_quadmesh_64times32_[0_1]_[1-2_1] elements vertices
% gS = polyMeshAuxStructure(vertices,elements);
% plotPolyMsh(gS)

% Darcy domain
left_D=0;right_D=1;bottom_D=0;top_D=1/2;
[M_D,T_D]=generate_quad_M_T(left_D,right_D,bottom_D,top_D,h,1);
elements = T_D';
vertices = M_D';
save D_quadmesh_64times32_[0_1]_[0_1-2] elements vertices
% gD = polyMeshAuxStructure(vertices,elements);
% plotPolyMsh(gD)

end % 
