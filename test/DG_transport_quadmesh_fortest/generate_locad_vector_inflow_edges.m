function r=generate_locad_vector_inflow_edges(func_name,t,DGT,DGE,DG_edge_flag,inflow_edges,basis_type,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vector_u_1)

r=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,1);

for i=1:length(inflow_edges)
    n=inflow_edges(i); % current edge
    ele=DGE(5,n);
    begin_point=DGE(1:2,n);end_point=DGE(3:4,n);

    temp=vector_u_1(1)*DG_edge_flag(2,n)+vector_u_1(2)*DG_edge_flag(3,n);
    
    b=zeros((basis_type+1)*(basis_type+2)/2,1);
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)*feval(func_name,begin_point(1),Gauss_point_local_1D(k),t)...
                *quadrilateral_local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,0,0)';
        end
    else  % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)*feval(func_name,Gauss_point_local_1D(k),begin_point(2),t)...
                *quadrilateral_local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,0,0)';
        end
    end
    
    for beta=1:(basis_type+1)*(basis_type+2)/2
        r((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,1)=r((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,1)+temp*b(beta,1);
    end
end

