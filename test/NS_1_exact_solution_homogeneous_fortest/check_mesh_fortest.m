function check_mesh_fortest
% 05/05/2016
clc
clearvars; close all;

format long e
Gauss_number_triangle=19;
Gauss_number_1D=10;
basis_type=2;
eipsilon=-1;  penalty=18;   co_beta=1;
% penalty=1, NIPG
% penalty=10, SIPG


Mesh_refinement=0;
[pressure_L2_error_1,speed_L2_error_1,speed_H1_error_1]=DG_Navier_Stokes_triangular_mesh(Mesh_refinement,... 
    Gauss_number_triangle,Gauss_number_1D,basis_type,eipsilon,penalty,co_beta);                                                                          
disp('complete_1');

time=toc

