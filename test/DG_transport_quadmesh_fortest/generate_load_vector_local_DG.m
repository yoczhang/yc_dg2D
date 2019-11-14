function r=generate_load_vector_local_DG(func_right_term,t,DGM,DGT...
                                        ,Gauss_coefficient_reference,Gauss_point_reference,basis_type,flag)

r=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,1);

for n=1:size(DGT,2)
    left_lower_point=DGM(:,DGT(1,n));
    h=[DGM(1,DGT(2,n))-DGM(1,DGT(1,n)) DGM(2,DGT(4,n))-DGM(2,DGT(1,n))];
    vertices=DGM(:,DGT(:,n));
    [Gauss_coefficient_local,Gauss_point_local]=generate_Gauss_local(Gauss_coefficient_reference,Gauss_point_reference,left_lower_point,h,vertices);
    
    b=zeros((basis_type+1)*(basis_type+2)/2,1);
    flag_ele=flag(1,n);
    for k=1:length(Gauss_coefficient_local)  
        a=feval(func_right_term,Gauss_point_local(k,1),Gauss_point_local(k,2),t,flag_ele)...
            *quadrilateral_local_basis(Gauss_point_local(k,1),Gauss_point_local(k,2),basis_type...
                                ,0,0);
        b=b+Gauss_coefficient_local(k)*a';
    end
    
    for i=1:(basis_type+1)*(basis_type+2)/2
        r((n-1)*(basis_type+1)*(basis_type+2)/2+i,1)=r((n-1)*(basis_type+1)*(basis_type+2)/2+i,1)+b(i,1);
    end
end