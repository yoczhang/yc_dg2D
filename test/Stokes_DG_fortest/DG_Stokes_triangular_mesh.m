function [pressure_L2_error,speed_L2_error,speed_H1_error]=DG_Stokes_triangular_mesh(Mesh_refinement,Gauss_number_triangle,Gauss_number_1D,basis_type,eipsilon,penalty,co_beta)
% 05/05/2016

nu=function_nu(0,0);
%% Triangulation
% Initial mesh
% DGM=[0,0;1,0;0,1;1,1;0,1;1,0]';
% DGT=[1,2,3;4,5,6]';


% [DGE,DG_edge_flag]=generate_DGE_DG_edge_flag(DGM,DGT);

%----------- yc test -------------------------
 %-------------------- Tri mesh ---------------------
left=0; right=1; bottom=1; top=2;
h_x = 1/2^(1+1); 
h_y = h_x;
h_partition=[h_x,h_y];
[node, elem] = generate_Tri_P_T(left,right,bottom,top,h_partition);
%-----------------------------------------------------

[DGM,DGT]=generate_DGM_DGT_DGE_DG_edge_flag(node',elem');
[DGE,DG_edge_flag]=generate_DGE_DG_edge_flag(DGM,DGT);
%-----------------------------------------------

% Begin uniform refinement
for i=1:Mesh_refinement
    [DGM,DGT]=uniformrefine_triangle(DGM,DGT);
end

[un_used,interior_edges]=find(DGE(6,:));
boundary_edges=setdiff(1:size(DGE,2),interior_edges);

[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(Gauss_number_triangle);
[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(Gauss_number_1D);

%% Assemble stiffness matrix

basis_speed=basis_type;
basis_pressure=basis_type-1;
A1=generate_stiffness_matrix_local_DG(DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_speed,basis_speed...
                                             ,1,0,1,0);
A2=generate_stiffness_matrix_local_DG(DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_speed,basis_speed...
                                             ,0,1,0,1);
[A3,A3_interE,A3_boundaryE]=generate_stiffness_matrix_local_DG_edge('integrand_1_on_edges','integrand_2_on_edges',DGM,DGT,DGE,DG_edge_flag...
                                                  ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                  ,basis_speed,eipsilon,penalty,co_beta,interior_edges,boundary_edges);
B1=generate_stiffness_matrix_local_DG(DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_pressure,basis_speed...
                                             ,0,0,1,0);
B2=generate_stiffness_matrix_local_DG(DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_pressure,basis_speed...
                                             ,0,0,0,1);
[B3,B4]=generate_stiffness_matrix_pressure_term_local_DG_edge(DGM,DGT,DGE,DG_edge_flag...
                                                                ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                                ,basis_pressure,basis_speed,interior_edges,boundary_edges);
%% Assemble load vector
b1=generate_load_vector_local_DG('function_f1',DGM,DGT...
                                        ,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_speed);
b2=generate_load_vector_local_DG('function_f2',DGM,DGT...
                                        ,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_speed);
                                    
b3=generate_load_vector_local_DG_boundary_edge('integrand_3_on_edges',DGM,DGT,DGE,DG_edge_flag...
                                                      ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                      ,basis_speed,boundary_edges,eipsilon,penalty,co_beta);
b4=generate_load_vector_local_DG_boundary_edge('integrand_4_on_edges',DGM,DGT,DGE,DG_edge_flag...
                                                      ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                      ,basis_speed,boundary_edges,eipsilon,penalty,co_beta);
b5=generate_load_vector_boundary_edge_divergence_term('function_u1','function_u2',DGT,DGE,DG_edge_flag,basis_pressure,boundary_edges...
                                                             ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D);                                                  

%% Solution
un_used_1=sparse(size(DGT,2)*(basis_speed+1)*(basis_speed+2)/2,size(DGT,2)*(basis_speed+1)*(basis_speed+2)/2);
% un_used_2=sparse(size(DGT,2)*(basis_pressure+1)*(basis_pressure+2)/2,size(DGT,2)*(basis_pressure+1)*(basis_pressure+2)/2);
% un_used_2=1e-6*eye(size(DGT,2)*(basis_pressure+1)*(basis_pressure+2)/2);
un_used_2=1e-6*generate_stiffness_matrix_local_DG(DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_pressure,basis_pressure...
                                             ,0,0,0,0);

left_matrix=[nu*(A1+A2+A3) un_used_1      -B1+B3;
             un_used_1     nu*(A1+A2+A3)  -B2+B4;
            (B1-B3)'      (B2-B4)'         un_used_2];
right_term=[b1+nu*b3;
            b2+nu*b4;
            -b5];
        
% %------- yc for test-----------------
% refA1 = A1; refA2 = A2; refA3 = A3; refA3_interE=A3_interE; refA3_boundaryE=A3_boundaryE;
% refB1 = B1; refB2 = B2; refB3 = B3; refB4 = B4; 
% refb1 = b1; refb2 = b2; refb3 = b3; refb4 = b4; refb5 = b5;
% ref_left_matrix = left_matrix; ref_right_term = right_term;
% save refStokesMatVec refA1 refA2 refA3 refB1 refB2 refB3 refB4 refb1 refb2 refb3 refb4 refb5 ...
%     ref_left_matrix ref_right_term refA3_interE refA3_boundaryE
% %--------------------------------------

% spy(left_matrix)
r=left_matrix\right_term;
clear A1 A2 A3 B1 B2 b1 b2 b3 b4 

%% Error
u1_h=r(1:size(DGT,2)*(basis_speed+1)*(basis_speed+2)/2,1);
u2_h=r(size(DGT,2)*(basis_speed+1)*(basis_speed+2)/2+1:2*size(DGT,2)*(basis_speed+1)*(basis_speed+2)/2,1);
p_h=r(2*size(DGT,2)*(basis_speed+1)*(basis_speed+2)/2+1:2*size(DGT,2)*(basis_speed+1)*(basis_speed+2)/2+size(DGT,2)*(basis_pressure+1)*(basis_pressure+2)/2,1);

pressure_L2_error=L2_H1_error(DGM,DGT,p_h,'function_p',basis_pressure,0,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle)

u1_L2_error=L2_H1_error(DGM,DGT,u1_h,'function_u1',basis_speed,0,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
u2_L2_error=L2_H1_error(DGM,DGT,u2_h,'function_u2',basis_speed,0,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
speed_L2_error=sqrt(u1_L2_error^2+u2_L2_error^2)

temp_1=L2_H1_error(DGM,DGT,u1_h,'function_u1_x',basis_speed,1,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
temp_2=L2_H1_error(DGM,DGT,u2_h,'function_u2_x',basis_speed,1,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
temp_3=L2_H1_error(DGM,DGT,u1_h,'function_u1_y',basis_speed,0,1,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
temp_4=L2_H1_error(DGM,DGT,u2_h,'function_u2_y',basis_speed,0,1,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
speed_H1_error=sqrt(temp_1^2+temp_2^2+temp_3^2+temp_4^2)

%% Figure
% uh=zeros(3*size(DGT,2),1);
% for n=1:size(DGT,2)
%     uh_local=p_h((n-1)*(basis_pressure+1)*(basis_pressure+2)/2+1:n*(basis_pressure+1)*(basis_pressure+2)/2,1);
%     uh(3*n-2,1)=fe_solution(DGM(1,DGT(1,n)),DGM(2,DGT(1,n)),uh_local,0,0,basis_pressure);
%     uh(3*n-1,1)=fe_solution(DGM(1,DGT(2,n)),DGM(2,DGT(2,n)),uh_local,0,0,basis_pressure);
%     uh(3*n,1)=fe_solution(DGM(1,DGT(3,n)),DGM(2,DGT(3,n)),uh_local,0,0,basis_pressure);
% end
% 
% figure(1);
% trisurf(DGT',DGM(1,:)',DGM(2,:)',uh);
% shading interp
% xlabel('x'), ylabel('y')
% colorbar 
% axis equal, axis off
% view(2)







