function temp=L2_H1_error(DGM,DGT,X,func_name,basis_type,der_order_x,der_order_y,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle)

temp=0;
for n=1:size(DGT,2)
    uh_local=X((n-1)*(basis_type+1)*(basis_type+2)/2+1:n*(basis_type+1)*(basis_type+2)/2,1);
    vertices=DGM(:,DGT(:,n));
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle...
                                                                                               ,Gauss_point_reference_triangle,vertices);
    for k=1:length(Gauss_coefficient_local_triangle)
        temp=temp+Gauss_coefficient_local_triangle(k)*(feval(func_name,Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2))...
           -fe_solution(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),uh_local,der_order_x,der_order_y,basis_type))^2;
    end
end
temp=sqrt(temp);