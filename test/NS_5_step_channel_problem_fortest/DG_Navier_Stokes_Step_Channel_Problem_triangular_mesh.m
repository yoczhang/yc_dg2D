function DG_Navier_Stokes_Step_Channel_Problem_triangular_mesh
% 05/05/2016

tic
%%
Gauss_number_triangle=19;Gauss_number_1D=10;
basis_type=2;eipsilon=1;penalty=1;co_beta=1;

%%
nu=function_nu(0,0);
basis_speed=basis_type;basis_pressure=basis_type-1;
N_loc_speed=(basis_speed+1)*(basis_speed+2)/2;N_loc_pressure=(basis_pressure+1)*(basis_pressure+2)/2;

%% Triangulation from FreeFem++
file=fopen('Stokes_nodes.dat');points=fscanf(file,'%f');fclose(file);
nppp=size(points);npp=nppp(1)/3;P=zeros(2,npp);
for i=1:npp
    P(1,i)=points(3*i-1);P(2,i)=points(3*i);
end
file=fopen('Stokes_elem.dat');connecticy=fscanf(file,'%f');fclose(file);
neee=size(connecticy);nee=neee(1)/9;T=zeros(3,nee);
for i=1:nee
    T(1,i)=connecticy(9*i-6)+1;T(2,i)=connecticy(9*i-3)+1;T(3,i)=connecticy(9*i)+1;
end

DGM=zeros(2,3*size(T,2));DGT=zeros(3,size(T,2));
for i=1:size(T,2)
    DGT(:,i)=[3*i-2 3*i-1 3*i]';
end

for i=1:size(T,2)
    for j=1:3
        DGM(:,DGT(j,i))=P(:,T(j,i));
    end
end

[DGE,DG_edge_flag]=generate_DGE_DG_edge_flag(DGM,DGT);

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
un_used_2=1e-10*M;
un_used_3=sparse(size(DGT,2)*N_loc_speed,size(DGT,2)*N_loc_pressure);
un_used_4=sparse(size(DGT,2)*N_loc_pressure,size(DGT,2)*N_loc_pressure);

left_matrix_0=[nu*(A1+A2+A3)   un_used_1        -B1+B3;
               un_used_1       nu*(A1+A2+A3)    -B2+B4;
              (B1-B3)'        (B2-B4)'           un_used_2];
right_term_0=[nu*b3;
              nu*b4;
             -b5];
r_0=left_matrix_0\right_term_0;

clear A1 A2 A3 B1 B2 B3 B4 b1 b2 b3 b4 b5

%% LOOP

err=1;
while(err>1e-8)
    u1_h=r_0(1:size(DGT,2)*N_loc_speed,1);
    u2_h=r_0(size(DGT,2)*N_loc_speed+1:2*size(DGT,2)*N_loc_speed,1);
    N1=generate_stiffness_matrix_nonlinear_term_1_2(DGM,DGT,basis_type,u1_h,u2_h,'integrand_nonlinear_function_1'...
                                              ,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    N2=generate_stiffness_matrix_nonlinear_term_3_4(DGT,DGE,DG_edge_flag,'integrand_nonlinear_function_2','integrand_nonlinear_function_3'...
                                                       ,'integrand_nonlinear_function_4','integrand_nonlinear_function_5',basis_type,u1_h,u2_h...
                                                       ,interior_edges,boundary_edges,Gauss_coefficient_reference_1D,Gauss_point_reference_1D);
    [b_N_1,b_N_2]=generate_load_vector_nonlinear_term_1_2(DGT,DGE,DG_edge_flag,boundary_edges,basis_type,'function_u1','function_u2'...
        ,u1_h,u2_h,Gauss_coefficient_reference_1D,Gauss_point_reference_1D);
    
    %---- yctest
    ref_N1 = N1;
%     ref_N2_interE = N2_interE;
%     ref_N2_boundaryE = N2_boundaryE;
    ref_b_N_1 = b_N_1;
    ref_b_N_2 = b_N_2;
    save nolinearTerm ref_N1 ref_b_N_1 ref_b_N_2
    %------------
    
                                                     
                                                     
    TEMP=[N1+N2 un_used_1 un_used_3;un_used_1 N1+N2 un_used_3;un_used_3' un_used_3' un_used_4];
    left_matrix=left_matrix_0+TEMP;
    right_term=right_term_0+[b_N_1;b_N_2;sparse(size(DGT,2)*N_loc_pressure,1)];
    r=left_matrix\right_term;
    
    err_1=generate_tolerate_error(DGM,DGT,u1_h,r(1:size(DGT,2)*N_loc_speed,1),basis_type,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    err_2=generate_tolerate_error(DGM,DGT,u2_h,r(size(DGT,2)*N_loc_speed+1:2*size(DGT,2)*N_loc_speed,1),basis_type,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle);
    err=sqrt(err_1^2+err_2^2)
    
    % Evaluate
    r_0=r;
end
clear N1 N2 u1_h u2_h b_N_1 b_N_2 TEMP left_matrix left_matrix_0 right_term right_term_0 r err err_1 err_2

    
%% Figure

uh_1=zeros(3*size(DGT,2),1);uh_2=zeros(3*size(DGT,2),1);ph=zeros(3*size(DGT,2),1);
for n=1:size(DGT,2) 
    uh_local_1=r_0((n-1)*N_loc_speed+1:n*N_loc_speed,1);
    uh_local_2=r_0(size(DGT,2)*N_loc_speed+(n-1)*N_loc_speed+1:size(DGT,2)*N_loc_speed+n*N_loc_speed,1);
    ph_local=r_0(2*size(DGT,2)*N_loc_speed+(n-1)*N_loc_pressure+1:2*size(DGT,2)*N_loc_speed+n*N_loc_pressure,1);

    uh_1(3*n-2,1)=fe_solution(DGM(1,DGT(1,n)),DGM(2,DGT(1,n)),uh_local_1,0,0,basis_speed);
    uh_1(3*n-1,1)=fe_solution(DGM(1,DGT(2,n)),DGM(2,DGT(2,n)),uh_local_1,0,0,basis_speed);
    uh_1(3*n,1)=fe_solution(DGM(1,DGT(3,n)),DGM(2,DGT(3,n)),uh_local_1,0,0,basis_speed);
    
    uh_2(3*n-2,1)=fe_solution(DGM(1,DGT(1,n)),DGM(2,DGT(1,n)),uh_local_2,0,0,basis_speed);
    uh_2(3*n-1,1)=fe_solution(DGM(1,DGT(2,n)),DGM(2,DGT(2,n)),uh_local_2,0,0,basis_speed);
    uh_2(3*n,1)=fe_solution(DGM(1,DGT(3,n)),DGM(2,DGT(3,n)),uh_local_2,0,0,basis_speed);
    
    ph(3*n-2,1)=fe_solution(DGM(1,DGT(1,n)),DGM(2,DGT(1,n)),ph_local,0,0,basis_pressure);
    ph(3*n-1,1)=fe_solution(DGM(1,DGT(2,n)),DGM(2,DGT(2,n)),ph_local,0,0,basis_pressure);
    ph(3*n,1)=fe_solution(DGM(1,DGT(3,n)),DGM(2,DGT(3,n)),ph_local,0,0,basis_pressure);
end
ph_min=min(ph);
ph=ph+abs(ph_min)+2;
figure(1000001);
trisurf(DGT',DGM(1,:)',DGM(2,:)',ph);
shading interp, xlabel('x'), ylabel('y'), colorbar
axis equal, axis off
title('Pressure');

norm_speed=sqrt(uh_1.^2+uh_2.^2);
figure(1000002);
trisurf(DGT',DGM(1,:)',DGM(2,:)',norm_speed);
shading interp, xlabel('x'), ylabel('y'), colorbar, view(2)
axis equal, axis off
title('Norm_speed');

%% Tecplot

% speed
f1=fopen('Ex9_case1_velocity_SD.dat','w'); fprintf(f1,'TITLE = ""\n'); fprintf(f1,'VARIABLES = "X" "Y" "U" "U1" "U2" \n'); 
total_number_of_nodes=length(norm_speed);total_number_of_elements=size(DGT,2);
fprintf(f1,'ZONE N=%d,E=%d\n',total_number_of_nodes,total_number_of_elements); fprintf(f1,'DATAPACKING=POINT,ZONETYPE=FETRIANGLE\n');
for i=1:total_number_of_nodes
   fprintf(f1,'%f\t%f\t%f\t%f\t%f\r\n',DGM(1,i),DGM(2,i),norm_speed(i),uh_1(i),uh_2(i));  
end
fprintf('\n');
for n=1:total_number_of_elements
   fprintf(f1,'%d\t%d\t%d\r\n',DGT(1,n),DGT(2,n),DGT(3,n)); 
end
fclose(f1);

% pressure
f1=fopen('Ex9_case1_pressure_SD.dat','w'); fprintf(f1,'TITLE = ""\n'); fprintf(f1,'VARIABLES = "X" "Y" "U"  \n'); 
total_number_of_nodes=length(ph);total_number_of_elements=size(DGT,2);
fprintf(f1,'ZONE N=%d,E=%d\n',total_number_of_nodes,total_number_of_elements); fprintf(f1,'DATAPACKING=POINT,ZONETYPE=FETRIANGLE\n');
for i=1:total_number_of_nodes
   fprintf(f1,'%f\t%f\t%f\t\r\n',DGM(1,i),DGM(2,i),ph(i));  
end
fprintf('\n');
for n=1:total_number_of_elements
   fprintf(f1,'%d\t%d\t%d\r\n',DGT(1,n),DGT(2,n),DGT(3,n)); 
end
fclose(f1);

time=toc



