function temp=L2_H1_error(DGM,DGT,X,t,func_name,basis_type,der_order_x,der_order_y,Gauss_coefficient_reference,Gauss_point_reference)

temp=0;
for n=1:size(DGT,2)
    uh_local=X((n-1)*(basis_type+1)*(basis_type+2)/2+1:n*(basis_type+1)*(basis_type+2)/2,1);
    left_lower_point=DGM(:,DGT(1,n));
    h=[DGM(1,DGT(2,n))-DGM(1,DGT(1,n)) DGM(2,DGT(4,n))-DGM(2,DGT(1,n))];
    vertices=DGM(:,DGT(:,n));
    [Gauss_coefficient_local,Gauss_point_local]=generate_Gauss_local(Gauss_coefficient_reference,Gauss_point_reference,left_lower_point,h,vertices);
    
    for k=1:length(Gauss_coefficient_local)
        temp=temp+Gauss_coefficient_local(k)*(feval(func_name,Gauss_point_local(k,1),Gauss_point_local(k,2),t)...
           -fe_solution(Gauss_point_local(k,1),Gauss_point_local(k,2),uh_local,der_order_x,der_order_y,basis_type))^2;
    end
end
temp=sqrt(temp);