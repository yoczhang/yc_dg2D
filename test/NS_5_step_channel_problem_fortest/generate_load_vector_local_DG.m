function r=generate_load_vector_local_DG(func_right_term,DGM,DGT...
                                        ,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type)

r=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,1);

for n=1:size(DGT,2)
    vertices=DGM(:,DGT(:,n));
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle...
                                                                                               ,Gauss_point_reference_triangle,vertices);
    b=zeros((basis_type+1)*(basis_type+2)/2,1);
    for k=1:length(Gauss_coefficient_local_triangle)  
        a=feval(func_right_term,Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2))...
            *local_basis(Gauss_point_local_triangle(k,1),Gauss_point_local_triangle(k,2),basis_type...
                                ,0,0);
        b=b+Gauss_coefficient_local_triangle(k)*a';
    end
    
    for i=1:(basis_type+1)*(basis_type+2)/2
        r((n-1)*(basis_type+1)*(basis_type+2)/2+i,1)=r((n-1)*(basis_type+1)*(basis_type+2)/2+i,1)+b(i,1);
    end
end