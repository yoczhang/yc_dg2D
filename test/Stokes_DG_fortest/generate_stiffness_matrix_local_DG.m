function r=generate_stiffness_matrix_local_DG(DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type_trial,basis_type_test...
                                             ,der_order_x_trial,der_order_y_trial,der_order_x_test,der_order_y_test)


r=sparse(size(DGT,2)*(basis_type_test+1)*(basis_type_test+2)/2,size(DGT,2)*(basis_type_trial+1)*(basis_type_trial+2)/2);

for n=1:size(DGT,2)
    vertices=DGM(:,DGT(:,n));
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle...
                                                                                               ,Gauss_point_reference_triangle,vertices);
    b=zeros((basis_type_trial+1)*(basis_type_trial+2)/2,(basis_type_test+1)*(basis_type_test+2)/2);
    for k=1:length(Gauss_coefficient_local_triangle)
        temp_trial=local_basis(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),basis_type_trial...
                                ,der_order_x_trial,der_order_y_trial);
        temp_test=local_basis(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),basis_type_test...
                                ,der_order_x_test,der_order_y_test);
        b=b+Gauss_coefficient_local_triangle(k)*temp_trial'*temp_test;
    end
    for i=1:(basis_type_trial+1)*(basis_type_trial+2)/2
        for j=1:(basis_type_test+1)*(basis_type_test+2)/2
            r((n-1)*(basis_type_test+1)*(basis_type_test+2)/2+j,(n-1)*(basis_type_trial+1)*(basis_type_trial+2)/2+i)...
                =r((n-1)*(basis_type_test+1)*(basis_type_test+2)/2+j,(n-1)*(basis_type_trial+1)*(basis_type_trial+2)/2+i)+b(i,j);
        end
    end
end


        