function check_polynomials

tic

Mesh_refinement=0;
eipsilon=-1;sigma=sqrt(2);beta=1;
delta_t=0.1;
theta=1; % Backward Euler


basis_type=1;
[L2_error_1,H1_error_1]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,16,8,theta);
disp('complete_1')

basis_type=2;
[L2_error_2,H1_error_2]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,16,8,theta);
disp('complete_2')

basis_type=3;
[L2_error_3,H1_error_3]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,16,8,theta);
disp('complete_3')

basis_type=4;
[L2_error_4,H1_error_4]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,16,8,theta);
disp('complete_4')

basis_type=5;
[L2_error_5,H1_error_5]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,100,10,theta);
disp('complete_5')

basis_type=6;
[L2_error_6,H1_error_6]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,100,10,theta);
disp('complete_6')

basis_type=7;
[L2_error_7,H1_error_7]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,100,10,theta);
disp('complete_7')

basis_type=8;
[L2_error_8,H1_error_8]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,100,10,theta);
disp('complete_8')

basis_type=9;
[L2_error_9,H1_error_9]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,100,10,theta);
disp('complete_9')

basis_type=10;
[L2_error_10,H1_error_10]=DG_transport_quadrilateral(Mesh_refinement,basis_type,eipsilon,sigma,beta,delta_t,100,10,theta);
disp('complete_10')

L2_error=[L2_error_1 L2_error_2 L2_error_3 L2_error_4 L2_error_5 L2_error_6 L2_error_7 L2_error_8 L2_error_9 L2_error_10]
H1_error=[H1_error_1 H1_error_2 H1_error_3 H1_error_4 H1_error_5 H1_error_6 H1_error_7 H1_error_8 H1_error_9 H1_error_10]

%% Figure
figure(1020);
N=zeros(1,10);
for i=1:10
    N(1,i)=(i+1)*(i+2);
end

semilogy(N,L2_error,'b-s',N,H1_error,'g-O','LineWidth',2);
% ,N,N.^(-5/2),'k--');
h1=legend('$||u-u_h||_0$','$||u-u_h||_{\cal E}$','Location','northeast');
set(h1,'Interpreter','latex')
xlabel('Number of degrees of freedom');
ylabel('Relative errors in various norms');

time=toc
