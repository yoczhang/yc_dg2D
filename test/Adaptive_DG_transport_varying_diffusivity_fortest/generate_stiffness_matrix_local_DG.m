function [r1,r2,r3]=generate_stiffness_matrix_local_DG(DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,eipsilon,basis_type,N_loc)

% r1 is mass matrix
% r2 is stiffness matrix
% r3 is matrix from convection term

r1=sparse(size(DGT,2)*N_loc,size(DGT,2)*N_loc);
r2=sparse(size(DGT,2)*N_loc,size(DGT,2)*N_loc);
r3=sparse(size(DGT,2)*N_loc,size(DGT,2)*N_loc);

for n=1:size(DGT,2)
    vertices=DGM(:,DGT(:,n));
    
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
    b1=zeros(N_loc,N_loc);
    b2=zeros(N_loc,N_loc);
    b3=zeros(N_loc,N_loc);
    for k=1:length(Gauss_coefficient_local_triangle)
        a1=triangular_local_basis(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),basis_type,0,0);
        a2=triangular_local_basis(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),basis_type,1,0);
        a3=triangular_local_basis(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),basis_type,0,1);
        
        b1=b1+Gauss_coefficient_local_triangle(k)*a1'*a1;
        b2=b2+eipsilon(1,n)*Gauss_coefficient_local_triangle(k)*(a2'*a2+a3'*a3);
        b3=b3+Gauss_coefficient_local_triangle(k)*a1'*a2;
    end
    for i=1:N_loc
        for j=1:N_loc
            r1((n-1)*N_loc+j,(n-1)*N_loc+i)=r1((n-1)*N_loc+j,(n-1)*N_loc+i)+b1(i,j);
            r2((n-1)*N_loc+j,(n-1)*N_loc+i)=r2((n-1)*N_loc+j,(n-1)*N_loc+i)+b2(i,j);
            r3((n-1)*N_loc+j,(n-1)*N_loc+i)=r3((n-1)*N_loc+j,(n-1)*N_loc+i)+b3(i,j);
        end
    end
end