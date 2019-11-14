function r=generate_stiffness_matrix_nonlinear_term_1_2(DGM,DGT,basis_type,u1_h,u2_h,func_name,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle)

r=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,size(DGT,2)*(basis_type+1)*(basis_type+2)/2);

for n=1:size(DGT,2)
    vertices=DGM(:,DGT(:,n));
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle...
                                                                                               ,Gauss_point_reference_triangle,vertices);
    b=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    solution_1=u1_h((n-1)*(basis_type+1)*(basis_type+2)/2+1:n*(basis_type+1)*(basis_type+2)/2,1);
    solution_2=u2_h((n-1)*(basis_type+1)*(basis_type+2)/2+1:n*(basis_type+1)*(basis_type+2)/2,1);
    
    for k=1:length(Gauss_coefficient_local_triangle)
        temp=feval(func_name,Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),solution_1,solution_2,basis_type);
        b=b+Gauss_coefficient_local_triangle(k)*temp;
    end
    
    for i=1:(basis_type+1)*(basis_type+2)/2
        for j=1:(basis_type+1)*(basis_type+2)/2
            r((n-1)*(basis_type+1)*(basis_type+2)/2+j,(n-1)*(basis_type+1)*(basis_type+2)/2+i)...
                =r((n-1)*(basis_type+1)*(basis_type+2)/2+j,(n-1)*(basis_type+1)*(basis_type+2)/2+i)+b(i,j);
        end
    end
end