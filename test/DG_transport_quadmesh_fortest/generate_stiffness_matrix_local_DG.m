function [r1,r2,r3]=generate_stiffness_matrix_local_DG...
    (DGM,DGT,Gauss_coefficient_reference,Gauss_point_reference,basis_type,flag,vector_u_1,vector_u_2)

% r1 is mass matrix
% r2 is stiffness matrix
% r3 is matrix from convection term

r1=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,size(DGT,2)*(basis_type+1)*(basis_type+2)/2);
r2=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,size(DGT,2)*(basis_type+1)*(basis_type+2)/2);
r3=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,size(DGT,2)*(basis_type+1)*(basis_type+2)/2);

for n=1:size(DGT,2)
    left_lower_point=DGM(:,DGT(1,n));
    h=[DGM(1,DGT(2,n))-DGM(1,DGT(1,n)) DGM(2,DGT(4,n))-DGM(2,DGT(1,n))];
    vertices=DGM(:,DGT(:,n));
    [Gauss_coefficient_local,Gauss_point_local]=generate_Gauss_local(Gauss_coefficient_reference,Gauss_point_reference,left_lower_point,h,vertices);
    b1=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    b2=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    b3=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    for k=1:length(Gauss_coefficient_local)
        a1=quadrilateral_local_basis(Gauss_point_local(k,1),Gauss_point_local(k,2),basis_type,0,0);
        a2=quadrilateral_local_basis(Gauss_point_local(k,1),Gauss_point_local(k,2),basis_type,1,0);
        a3=quadrilateral_local_basis(Gauss_point_local(k,1),Gauss_point_local(k,2),basis_type,0,1);
        
        b1=b1+Gauss_coefficient_local(k)*a1'*a1;
        b2=b2+Gauss_coefficient_local(k)*a2'*a2+Gauss_coefficient_local(k)*a3'*a3;
        if flag(n)==1
            b3=b3+vector_u_1(1)*Gauss_coefficient_local(k)*a1'*a2+vector_u_1(2)*Gauss_coefficient_local(k)*a1'*a3;
        else
            b3=b3+vector_u_2(1)*Gauss_coefficient_local(k)*a1'*a2+vector_u_2(2)*Gauss_coefficient_local(k)*a1'*a3;
        end
    end
    for i=1:(basis_type+1)*(basis_type+2)/2
        for j=1:(basis_type+1)*(basis_type+2)/2
            r1((n-1)*(basis_type+1)*(basis_type+2)/2+j,(n-1)*(basis_type+1)*(basis_type+2)/2+i)...
                =r1((n-1)*(basis_type+1)*(basis_type+2)/2+j,(n-1)*(basis_type+1)*(basis_type+2)/2+i)+b1(i,j);
            r2((n-1)*(basis_type+1)*(basis_type+2)/2+j,(n-1)*(basis_type+1)*(basis_type+2)/2+i)...
                =r2((n-1)*(basis_type+1)*(basis_type+2)/2+j,(n-1)*(basis_type+1)*(basis_type+2)/2+i)+b2(i,j);
            r3((n-1)*(basis_type+1)*(basis_type+2)/2+j,(n-1)*(basis_type+1)*(basis_type+2)/2+i)...
                =r3((n-1)*(basis_type+1)*(basis_type+2)/2+j,(n-1)*(basis_type+1)*(basis_type+2)/2+i)+b3(i,j);
        end
    end
end


        