function r=generate_stiffness_matrix_from_convection_term_2(DGE,DGT,DG_edge_flag,basis_type,outflow_edges,flag,vector_u_1,vector_u_2,integrand_func_name...
                                                           ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D)

% Because we know the orientation of the velocity (-0.2,-0.1), the upwind value of
% concentration must be taken as the limitation of that on the larger element

r=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,size(DGT,2)*(basis_type+1)*(basis_type+2)/2);

%% 
for j=1:length(outflow_edges)
    i=outflow_edges(j);
    ele=DGE(5,i);
    begin_point=DGE(1:2,i);end_point=DGE(3:4,i);
    
    b=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    if flag(ele)==1
        vector_u=vector_u_1;
    else
        vector_u=vector_u_2;
    end
    temp=vector_u(1)*DG_edge_flag(2,i)+vector_u(2)*DG_edge_flag(3,i);
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+temp*Gauss_coefficient_local_1D(k)*feval(integrand_func_name,begin_point(1),Gauss_point_local_1D(k),basis_type);
        end
    else % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+temp*Gauss_coefficient_local_1D(k)*feval(integrand_func_name,Gauss_point_local_1D(k),begin_point(2),basis_type);
        end
    end
    
    for alpha=1:(basis_type+1)*(basis_type+2)/2
        for beta=1:(basis_type+1)*(basis_type+2)/2
            r((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele-1)*(basis_type+1)*(basis_type+2)/2+alpha)...
                =r((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele-1)*(basis_type+1)*(basis_type+2)/2+alpha)+b(alpha,beta);
        end
    end
end

