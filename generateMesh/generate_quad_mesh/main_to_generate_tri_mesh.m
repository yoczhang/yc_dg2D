function main_to_generate_tri_mesh
%
%   Here we make the rules for naming the mesh, 
%   such as: 'S_quadmesh_4times2_[0_1]_[1-2_1]',
%       S_quadmesh, is the mesh name,
%       4time2, is the x-coordinate and y-coordinate elem arrange,
%       [0_1], is the x-coordinate domain [0,1]
%       [1-2_1], is the y-coordinate domain [1/2,1], and 1-2 is the transformation of 1/2.
%       
%
%   YcZhang 24/9/2017
%
%   Last modified 24/9/2017
%

clear
close all

% generate the single Tri mesh
getPT_TriMesh;


end % function


%% ---- sub function 1
function getPT_TriMesh

left=0; right=2; bottom=0; top=1;
h_x = 1/16; h_y = h_x;
h=[h_x,h_y];

[node, elem] = generate_Tri_P_T(left,right,bottom,top,h);
% becasuse 
save Trimesh_2times32times16_[0_2]_[0_1] node elem
end % function
