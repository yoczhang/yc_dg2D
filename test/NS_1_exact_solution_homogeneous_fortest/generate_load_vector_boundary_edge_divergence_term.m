function r=generate_load_vector_boundary_edge_divergence_term(func_u1,func_u2,DGT,DGE,DG_edge_flag,basis_type_test,boundary_edges...
                                                             ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D)

r=sparse(size(DGT,2)*(basis_type_test+1)*(basis_type_test+2)/2,1);

%% Boundary edges
for i=1:length(boundary_edges)
    n=boundary_edges(i);ele=DGE(5,n);
    begin_point=DGE(1:2,n);end_point=DGE(3:4,n);
    
    b=zeros(1,(basis_type_test+1)*(basis_type_test+2)/2);
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)...
                     *local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type_test,0,0)...
                     *(feval(func_u1,begin_point(1),Gauss_point_local_1D(k))*DG_edge_flag(2,n)...
                      +feval(func_u2,begin_point(1),Gauss_point_local_1D(k))*DG_edge_flag(3,n));             
        end
    else  % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)...
                     *local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type_test,0,0)...
                     *(feval(func_u1,Gauss_point_local_1D(k),begin_point(2))*DG_edge_flag(2,n)...
                      +feval(func_u2,Gauss_point_local_1D(k),begin_point(2))*DG_edge_flag(3,n));   
        end
    end
    
    for beta=1:(basis_type_test+1)*(basis_type_test+2)/2
        r((ele-1)*(basis_type_test+1)*(basis_type_test+2)/2+beta,1)...
                =r((ele-1)*(basis_type_test+1)*(basis_type_test+2)/2+beta,1)...
                  +b(1,beta);
    end
end
