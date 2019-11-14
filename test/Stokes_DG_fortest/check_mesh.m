function check_mesh
% 05/05/2016
clc
close all
clearvars
tic
format long e
Gauss_number_triangle=19;
Gauss_number_1D=10;
basis_type=2;
eipsilon=-1;  penalty=18;   co_beta=1;


Mesh_refinement=0;
[pressure_L2_error_0,speed_L2_error_0,speed_H1_error_0]=DG_Stokes_triangular_mesh(Mesh_refinement,Gauss_number_triangle,Gauss_number_1D...
                                                                           ,basis_type,eipsilon,penalty,co_beta);                                                                       
disp('complete_0');
Mesh_refinement=1;
[pressure_L2_error_1,speed_L2_error_1,speed_H1_error_1]=DG_Stokes_triangular_mesh(Mesh_refinement,Gauss_number_triangle,Gauss_number_1D...
                                                                           ,basis_type,eipsilon,penalty,co_beta);                                                                       
disp('complete_1');
Mesh_refinement=2;
[pressure_L2_error_2,speed_L2_error_2,speed_H1_error_2]=DG_Stokes_triangular_mesh(Mesh_refinement,Gauss_number_triangle,Gauss_number_1D...
                                                                           ,basis_type,eipsilon,penalty,co_beta);                                                                       
disp('complete_2');
Mesh_refinement=3;
[pressure_L2_error_3,speed_L2_error_3,speed_H1_error_3]=DG_Stokes_triangular_mesh(Mesh_refinement,Gauss_number_triangle,Gauss_number_1D...
                                                                           ,basis_type,eipsilon,penalty,co_beta);          
disp('complete_3');

pressure_L2_error=[pressure_L2_error_0,pressure_L2_error_1,pressure_L2_error_2,pressure_L2_error_3]
speed_L2_error=[speed_L2_error_0,speed_L2_error_1,speed_L2_error_2,speed_L2_error_3]
speed_H1_error=[speed_H1_error_0,speed_H1_error_1,speed_H1_error_2,speed_H1_error_3]

pressure_L2_error_order=zeros(1,3);
speed_L2_error_order=zeros(1,3);
speed_H1_error_order=zeros(1,3);
for i=1:3
    pressure_L2_error_order(1,i)=log(pressure_L2_error(i)/pressure_L2_error(i+1))/log(2);
    speed_L2_error_order(1,i)=log(speed_L2_error(i)/speed_L2_error(i+1))/log(2);
    speed_H1_error_order(1,i)=log(speed_H1_error(i)/speed_H1_error(i+1))/log(2);
end
pressure_L2_error_order,speed_L2_error_order,speed_H1_error_order
time=toc

