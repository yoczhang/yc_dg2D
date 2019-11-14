function r=generate_stiffness_matrix_local_DG_edge(integrand_func_name_1,DGT,DGE,DG_edge_flag,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                  ,basis_type,eipsilon,penalty_1,co_beta,interior_edges)


r=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,size(DGT,2)*(basis_type+1)*(basis_type+2)/2);

%% Interior edges
for i=1:length(interior_edges)
    n=interior_edges(i);ele=DGE(5,n);ele_neighbor=DGE(5,DGE(6,n));
    begin_point=DGE(1:2,n);end_point=DGE(3:4,n);
    edge=sqrt((end_point(1)-begin_point(1)).^2+(end_point(2)-begin_point(2)).^2);
    
    b=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    b_2=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)*feval(integrand_func_name_1,begin_point(1),Gauss_point_local_1D(k)...
                ,penalty_1,DG_edge_flag(2:4,n),DG_edge_flag(2:4,n),basis_type,eipsilon,edge,co_beta);
            b_2=b_2+Gauss_coefficient_local_1D(k)*feval(integrand_func_name_1,begin_point(1),Gauss_point_local_1D(k)...
                ,penalty_1,DG_edge_flag(2:4,n),DG_edge_flag(2:4,DGE(6,n)),basis_type,eipsilon,edge,co_beta);
        end
    elseif begin_point(2)==end_point(2) % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)*feval(integrand_func_name_1,Gauss_point_local_1D(k),begin_point(2)...
                ,penalty_1,DG_edge_flag(2:4,n),DG_edge_flag(2:4,n),basis_type,eipsilon,edge,co_beta);
            b_2=b_2+Gauss_coefficient_local_1D(k)*feval(integrand_func_name_1,Gauss_point_local_1D(k),begin_point(2)...
                ,penalty_1,DG_edge_flag(2:4,n),DG_edge_flag(2:4,DGE(6,n)),basis_type,eipsilon,edge,co_beta);
        end
    else
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        slope=(end_point(2)-begin_point(2))/(end_point(1)-begin_point(1));
        Jacobi=sqrt(1+slope^2);
        for k=1:length(Gauss_coefficient_local_1D)
            x=Gauss_point_local_1D(k);
            y=slope*(x-begin_point(1))+begin_point(2);
            b=b+Gauss_coefficient_local_1D(k)*Jacobi*feval(integrand_func_name_1,x,y...
                ,penalty_1,DG_edge_flag(2:4,n),DG_edge_flag(2:4,n),basis_type,eipsilon,edge,co_beta);
            b_2=b_2+Gauss_coefficient_local_1D(k)*Jacobi*feval(integrand_func_name_1,x,y...
                ,penalty_1,DG_edge_flag(2:4,n),DG_edge_flag(2:4,DGE(6,n)),basis_type,eipsilon,edge,co_beta);
        end
    end
    
    for alpha=1:(basis_type+1)*(basis_type+2)/2
        for beta=1:(basis_type+1)*(basis_type+2)/2
            r((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele-1)*(basis_type+1)*(basis_type+2)/2+alpha)...
                =r((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele-1)*(basis_type+1)*(basis_type+2)/2+alpha)+b(alpha,beta);
            r((ele_neighbor-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele-1)*(basis_type+1)*(basis_type+2)/2+alpha)...
                =r((ele_neighbor-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele-1)*(basis_type+1)*(basis_type+2)/2+alpha)+b_2(alpha,beta);
        end
    end
end

end % function
