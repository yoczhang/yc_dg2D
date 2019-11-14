function [L2_error,H1_error]=DG_transport_quadrilateral...
    (Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,Gauss_point_number_1,Gauss_point_number_2,theta)
% 04/28/2016


%% Parameters
% D=1;phi=0.1;T=0.5;
D=1;phi=1;T=0.5;

%% Using uniform quadrilateral partition, initial triangulation
% DGM=[0,0;10,0;10,7;0,3;0,3;10,7;10,10;0,10]';
% DGT=[1,2,3,4;5,6,7,8]';
DGM = [0,0; 1,0; 1,1; 0,1]';
DGT = [1,2,3,4]';




%% Begin uniform refinement
Mesh_refinement = Mesh_refinement+1;
for i=1:Mesh_refinement
    [DGM,DGT]=uniformrefine_quadrilateral(DGM,DGT);
end

% %---------------------------------------------------
% load('quadmesh_16elem');
% [DGM,DGT]=generate_DGM_DGT(vertices',elements');
% % trisurf(DGT',DGM(1,:)',DGM(2,:)',zeros(size(DGM,2),1),'FaceColor','w','EdgeColor','b');
% % view(2);
% % g = polyMeshAuxStructure(DGM, DGT);
% % plotPolyMsh(g)
% %-------------------------------------------------------

[DGE,DG_edge_flag]=generate_DGE_DG_edge_flag(DGM,DGT);

[un_used,interior_edges]=find(DGE(6,:));
boundary_edges=setdiff(1:size(DGE,2),interior_edges);

% Inflow and outflow edges
% [un_used,temp1]=find(DG_edge_flag(2,boundary_edges(:))==1);
middle_point_x=(DGE(1,boundary_edges(:))+DGE(3,boundary_edges(:)))/2;
middle_point_y=(DGE(2,boundary_edges(:))+DGE(4,boundary_edges(:)))/2;
[un_used,temp2]=find(middle_point_x==1);
[un_used,temp3]=find(middle_point_y==1);
temp=[temp2,temp3];
inflow_edges=boundary_edges(temp(:));
outflow_edges=setdiff(boundary_edges,inflow_edges);
clear temp1 temp2 temp middle_point_y un_used

% Gaussian quadrature points and weights in reference element
[Gauss_coefficient_reference,Gauss_point_reference]=generate_Gauss_reference(Gauss_point_number_1);
[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(Gauss_point_number_2);

% Seperate the domain
flag=ones(1,size(DGT,2)); % 1 or 2 denote the 1st or 2nd subdomain

% Take different vector_u in different subdomain
vector_u_1=[-1 -0.4];vector_u_2=[0 0];
% vector_u_1=[1, 0.5];vector_u_2=[0 0];


%% Assemble stiffness matrix
[M,A1,A2]=generate_stiffness_matrix_local_DG...
    (DGM,DGT,Gauss_coefficient_reference,Gauss_point_reference,basis_type,flag,vector_u_1,vector_u_2);

A3=generate_stiffness_matrix_local_DG_edge('integrand_1_on_edges',DGT,DGE,DG_edge_flag,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                  ,basis_type,eipsilon,sigma,beta,interior_edges);
C1=generate_stiffness_matrix_from_convection_term_1(DGE,DGT,DG_edge_flag,basis_type,interior_edges,flag,vector_u_1,vector_u_2,'integrand_4_on_edges'...
                                                           ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D);
C2=generate_stiffness_matrix_from_convection_term_2(DGE,DGT,DG_edge_flag,basis_type,outflow_edges,flag,vector_u_1,vector_u_2,'integrand_4_on_edges'...
                                                           ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D);

%% LOOP
% Because the load vector depends on the time, we must write it into the time cycle

N_T=T/delta_t;

left_matrix=phi*M+delta_t*theta*(D*A1-A2+A3+C1+C2);
right_matrix=phi*M-delta_t*(1-theta)*(D*A1-A2+A3+C1+C2);

% The L2 projection of the initial analytical solution
b_0=generate_load_vector_local_DG_initial_solution('accurate_solution',0,DGM,DGT...
                                        ,Gauss_coefficient_reference,Gauss_point_reference,basis_type);
r_0=M\b_0;

% %------------------------------------
% M_ref = M;
% A1_ref = A1;
% A2_ref = A2;
% A3_ref = A3;
% C1_ref = C1;
% C2_ref = C2;
% left_matrix_ref = left_matrix;
% right_matrix_ref = right_matrix;
% b_0_ref = b_0;
% r_0_ref = r_0;
% save transport_mat_ref_16elem M_ref A1_ref A2_ref A3_ref C1_ref C2_ref left_matrix_ref right_matrix_ref b_0_ref r_0_ref
% %---------------------------------------

clear M A1 A2 A3 C1 C2 b_0


for n=0:N_T-1
    t_1=n*delta_t; % current time
    t_2=(n+1)*delta_t; % next time
    t=theta*t_2+(1-theta)*t_1;
    
    % Assemble the load vector
    b1=generate_load_vector_local_DG('function_f',t,DGM,DGT...
                                        ,Gauss_coefficient_reference,Gauss_point_reference,basis_type,flag);
    b2=generate_locad_vector_inflow_edges('function_inflow',t,DGT,DGE,DG_edge_flag,inflow_edges,basis_type...
                                        ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vector_u_1);
    
    % Solve
    r=left_matrix\(right_matrix*r_0+delta_t*(b1-b2));
    
%     %-------------------------------------------------------
%     b1_ref = b1; b2_ref = b2; r_ref = r; t_ref=t;
%     save transport_loop_ref_16elem b1_ref b2_ref r_ref t_ref
%     %------------------------------------------------------

    % Evaluate
    r_0=r;
    
end

%% Error
L2_error=L2_H1_error(DGM,DGT,r,T,'accurate_solution',basis_type,0,0,Gauss_coefficient_reference,Gauss_point_reference)/sqrt(125/4);

temp_1=L2_H1_error(DGM,DGT,r,T,'accurate_solution_x',basis_type,1,0,Gauss_coefficient_reference,Gauss_point_reference)/sqrt(pi^2/2);
temp_2=L2_H1_error(DGM,DGT,r,T,'accurate_solution_y',basis_type,0,1,Gauss_coefficient_reference,Gauss_point_reference)/sqrt(pi^2/2);
H1_error=sqrt(temp_1^2+temp_2^2);

                                             

                                    

