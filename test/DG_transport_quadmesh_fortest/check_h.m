function check_h
clearvars;
close all;
clc;
tic

basis_type=1;
if basis_type == 1
    sigma = 6;
elseif basis_type == 2
    sigma = 18;
end
eipsilon=-1;
% eipsilon=1;sigma=1;
beta=1;
delta_t=0.01;
Gauss_point_number_1=16;
Gauss_point_number_2=8;
theta=1; % Backward Euler


Mesh_refinement=0;
[L2_error_0,H1_error_0]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,Gauss_point_number_1,Gauss_point_number_2,theta);
disp('complete_0')

Mesh_refinement=1;
[L2_error_1,H1_error_1]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,Gauss_point_number_1,Gauss_point_number_2,theta);
disp('complete_1')

Mesh_refinement=2;
[L2_error_2,H1_error_2]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,Gauss_point_number_1,Gauss_point_number_2,theta);
disp('complete_2')

Mesh_refinement=3;
[L2_error_3,H1_error_3]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,Gauss_point_number_1,Gauss_point_number_2,theta);
disp('complete_3')

Mesh_refinement=4;
[L2_error_4,H1_error_4]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,Gauss_point_number_1,Gauss_point_number_2,theta);
disp('complete_4')

L2_error=[L2_error_0 L2_error_1 L2_error_2 L2_error_3 L2_error_4]
H1_error=[H1_error_0 H1_error_1 H1_error_2 H1_error_3 H1_error_4]

%% Figure
figure(1020);
N=[2, 2*4, 2*4^2, 2*4^3, 2*4^4];

loglog(N,L2_error,'b-s',N,H1_error,'g-O',N,N.^(-1.5),'r--',N,N.^(-1),'r-.','LineWidth',2);
% ,N,N.^(-5/2),'k--');
h1=legend('$||u-u_h||_0$','$||u-u_h||_{\cal E}$','slope=1.5','slope=1','Location','northeast');
set(h1,'Interpreter','latex')
xlabel('Number of degrees of freedom');
ylabel('Relative errors in various norms');




time=toc

