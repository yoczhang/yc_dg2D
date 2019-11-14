function temp=generate_tolerate_error(DGM,DGT,uh_new,uh_old,basis_type,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle)

temp=0;
for n=1:size(DGT,2)
    solution_1=uh_new((n-1)*(basis_type+1)*(basis_type+2)/2+1:n*(basis_type+1)*(basis_type+2)/2,1);
    solution_2=uh_old((n-1)*(basis_type+1)*(basis_type+2)/2+1:n*(basis_type+1)*(basis_type+2)/2,1);
    
    vertices=DGM(:,DGT(:,n));
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle...
                                                                                               ,Gauss_point_reference_triangle,vertices);
    for k=1:length(Gauss_coefficient_local_triangle)
        temp=temp+Gauss_coefficient_local_triangle(k)...
            *(fe_solution(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),solution_1,0,0,basis_type)...
           -fe_solution(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),solution_2,0,0,basis_type))^2;
    end
end
temp=sqrt(temp);