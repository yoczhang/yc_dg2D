function [pressure_L2_error, speed_L2_error, speed_H1_error]=DG_Navier_Stokes_triangular_mesh(Mesh_refinement,Gauss_number_triangle,Gauss_number_1D,basis_type,eipsilon,penalty,co_beta)
% 05/05/2016


nu=function_nu(0,0);
basis_speed=basis_type;
basis_pressure=basis_type-1;
disp_speed_name = ['basis_speed = ',num2str(basis_speed),'.'];
disp(disp_speed_name);
disp_pressure_name = ['basis_pressure = ',num2str(basis_pressure),'.'];
disp(disp_pressure_name);

N_loc_speed=(basis_speed+1)*(basis_speed+2)/2;N_loc_pressure=(basis_pressure+1)*(basis_pressure+2)/2;
%% Triangulation
% Initial mesh
% DGM=[0,0;1,0;0,1;1,1;0,1;1,0]';
% DGT=[1,2,3;4,5,6]';
% 
% % Begin uniform refinement
% for i=1:Mesh_refinement
%     [DGM,DGT]=uniformrefine_triangle(DGM,DGT);
% end
% 
% [DGE,DG_edge_flag]=generate_DGE_DG_edge_flag(DGM,DGT);

%----------------------- yc test ---------------------------

%-- mesh 1
% [node,elem] = squaremesh([0,2,0,1],1/20);
% bdFlag = setboundary(node,elem,'Dirichlet');
% for k = 1:1
% 	[node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
% end

%-- mesh 2
% h_x = 1/20; h_y = h_x;
% h=[h_x,h_y];
% left=0;right=2;bottom=0;top=1;
% [M,T]=generate_quad_P_T(left,right,bottom,top,h,1);
% node = M'; elem = T';

%-- mesh 3
[node, elem] = get_BR_paper_mesh(1/1);

%-- mesh 4
% [node, elem] = triMeshOnPolygonDomain([0,2,2,0],[0,0,1,1],1/20);
%load('Delaunaymesh_80times40_[0_2]_[0_1]');

[DGM,DGT]=generate_DGM_DGT_DGE_DG_edge_flag(node',elem');
[DGE,DG_edge_flag]=generate_DGE_DG_edge_flag(DGM,DGT);
%-------------------------------------------------------------

[un_used,interior_edges]=find(DGE(6,:));
boundary_edges=setdiff(1:size(DGE,2),interior_edges);

[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(Gauss_number_triangle);
[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(Gauss_number_1D);

%% Assemble stiffness matrix
M=generate_stiffness_matrix_local_DG(DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_pressure,basis_pressure...
                                             ,0,0,0,0);
A1=generate_stiffness_matrix_local_DG(DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_speed,basis_speed...
                                             ,1,0,1,0);
A2=generate_stiffness_matrix_local_DG(DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_speed,basis_speed...
                                             ,0,1,0,1);
A3=generate_stiffness_matrix_local_DG_edge('integrand_1_on_edges','integrand_2_on_edges',DGM,DGT,DGE,DG_edge_flag...
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

%% Initial solution---solution of Stokes equation
un_used_1=sparse(size(DGT,2)*N_loc_speed,size(DGT,2)*N_loc_speed);
% un_used_2=1e-6*eye(size(DGT,2)*N_loc_pressure);
un_used_2=1e-6*M;
un_used_3=sparse(size(DGT,2)*N_loc_speed,size(DGT,2)*N_loc_pressure);
un_used_4=sparse(size(DGT,2)*N_loc_pressure,size(DGT,2)*N_loc_pressure);

left_matrix_0=[nu*(A1+A2+A3)   un_used_1        -B1+B3;
               un_used_1       nu*(A1+A2+A3)    -B2+B4;
              (B1-B3)'        (B2-B4)'           un_used_2];
right_term_0=[b1+nu*b3;
              b2+nu*b4;
             -b5];
r_0=left_matrix_0\right_term_0;

ref_A1 = A1; ref_A2 = A2; ref_A3 = A3;
ref_B1 = B1; ref_B2 = B2; ref_B3 = B3; ref_B4 = B4;
ref_b1 = b1; ref_b2 = b2; ref_b3 = b3; ref_b4 = b4; ref_b5 = b5; 
ref_mat0 = left_matrix_0;
ref_right0 = right_term_0;

save ref_mat ref_A1 ref_A2 ref_A3 ref_B1 ref_B2 ref_B3 ref_B4 ref_b1 ref_b2 ref_b3 ref_b4 ref_b5 ref_mat0 ref_right0
clear A1 A2 A3 B1 B2 B3 B4 b1 b2 b3 b4 b5
clear ref_A1 ref_A2 ref_A3 ref_B1 ref_B2 ref_B3 ref_B4 ref_b1 ref_b2 ref_b3 ref_b4 ref_b5 ref_mat0 ref_right0

%% LOOP

err=1;
iterationStep = 0;
n1 = 0; % control the disp()
Nstep = 99;
while(err>1e-10 && iterationStep<=99)
    u1_h=r_0(1:size(DGT,2)*N_loc_speed,1);
    u2_h=r_0(size(DGT,2)*N_loc_speed+1:2*size(DGT,2)*N_loc_speed,1);
    N1=generate_stiffness_matrix_nonlinear_term_1_2(DGM,DGT,basis_type,u1_h,u2_h,'integrand_nonlinear_function_1'...
                                              ,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    [N2_interE, N2_boundaryE]=generate_stiffness_matrix_nonlinear_term_3_4(DGT,DGE,DG_edge_flag,'integrand_nonlinear_function_2','integrand_nonlinear_function_3'...
                                                       ,'integrand_nonlinear_function_4','integrand_nonlinear_function_5',basis_type,u1_h,u2_h...
                                                       ,interior_edges,boundary_edges,Gauss_coefficient_reference_1D,Gauss_point_reference_1D);
                                                   
    [b_N_1,b_N_2]=generate_load_vector_nonlinear_term_1_2(DGT,DGE,DG_edge_flag,boundary_edges,basis_type,'function_u1','function_u2',...
        u1_h,u2_h,Gauss_coefficient_reference_1D,Gauss_point_reference_1D);
    N2=N2_interE+N2_boundaryE;
    TEMP=[N1+N2 un_used_1 un_used_3;un_used_1 N1+N2 un_used_3;un_used_3' un_used_3' un_used_4];
    left_matrix=left_matrix_0+TEMP;
    right_term=right_term_0+[b_N_1;b_N_2;sparse(size(DGT,2)*N_loc_pressure,1)];
    
    %---- yctest
    ref_N1 = N1;
    ref_N2_interE = N2_interE;
    ref_N2_boundaryE = N2_boundaryE;
    ref_b_N_1 = b_N_1;
    ref_b_N_2 = b_N_2;
    save nolinearTerm ref_N1 ref_N2_interE ref_N2_boundaryE ref_b_N_1 ref_b_N_2
    %------------
    
    r=left_matrix\right_term;
    
    err_1=generate_tolerate_error(DGM,DGT,u1_h,r(1:size(DGT,2)*N_loc_speed,1),basis_type,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    err_2=generate_tolerate_error(DGM,DGT,u2_h,r(size(DGT,2)*N_loc_speed+1:2*size(DGT,2)*N_loc_speed,1),basis_type,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    err=sqrt(err_1^2+err_2^2);
    
    % Evaluate
    r_0=r;
    
    % iterationStep
    iterationStep = iterationStep + 1;
    % iterationStep
    iterationStep = iterationStep + 1;
    if iterationStep==0 || iterationStep>=n1
        dispname1 = ['in the ', num2str(iterationStep),'-th iteration.'];
        disp(dispname1)
        n1 = n1 + Nstep/10;
    end % if
end


    
%% Error
u1_h=r(1:size(DGT,2)*N_loc_speed,1);
u2_h=r(size(DGT,2)*N_loc_speed+1:2*size(DGT,2)*N_loc_speed,1);
p_h=r(2*size(DGT,2)*N_loc_speed+1:2*size(DGT,2)*N_loc_speed+size(DGT,2)*N_loc_pressure,1);

pressure_L2_error=L2_H1_error(DGM,DGT,p_h,'function_p',basis_pressure,0,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);

u1_L2_error=L2_H1_error(DGM,DGT,u1_h,'function_u1',basis_speed,0,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
u2_L2_error=L2_H1_error(DGM,DGT,u2_h,'function_u2',basis_speed,0,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
speed_L2_error=sqrt(u1_L2_error^2+u2_L2_error^2);

temp_1=L2_H1_error(DGM,DGT,u1_h,'function_u1_x',basis_speed,1,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
temp_2=L2_H1_error(DGM,DGT,u2_h,'function_u2_x',basis_speed,1,0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
temp_3=L2_H1_error(DGM,DGT,u1_h,'function_u1_y',basis_speed,0,1,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
temp_4=L2_H1_error(DGM,DGT,u2_h,'function_u2_y',basis_speed,0,1,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
speed_H1_error=sqrt(temp_1^2+temp_2^2+temp_3^2+temp_4^2);







