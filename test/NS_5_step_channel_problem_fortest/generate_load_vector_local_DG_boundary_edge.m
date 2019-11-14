function r=generate_load_vector_local_DG_boundary_edge(integrand_func_name,DGM,DGT,DGE,DG_edge_flag...
                                                      ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                      ,basis_type,boundary_edges,eipsilon,penalty,co_beta)

r=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,1);

for i=1:length(boundary_edges)
    n=boundary_edges(i);ele=DGE(5,n);
    begin_point=DGE(1:2,n);end_point=DGE(3:4,n);
    edge=sqrt((end_point(1)-begin_point(1)).^2+(end_point(2)-begin_point(2)).^2);
    
    b=zeros((basis_type+1)*(basis_type+2)/2,1);
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)*feval(integrand_func_name,begin_point(1),Gauss_point_local_1D(k)...
                ,penalty,DG_edge_flag(2:4,n),basis_type,eipsilon,edge,co_beta);
        end
    else  % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)*feval(integrand_func_name,Gauss_point_local_1D(k),begin_point(2)...
                ,penalty,DG_edge_flag(2:4,n),basis_type,eipsilon,edge,co_beta);
        end
    end
    
    for beta=1:(basis_type+1)*(basis_type+2)/2
        r((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,1)=r((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,1)+b(beta,1);
    end
end