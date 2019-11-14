function r=generate_stiffness_matrix_diffusion_upwind_scheme_alpha_part(integrand_func_name_1,DGT,DGE,DG_edge_flag,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                  ,F_HP_i,basis_type,N_loc,eipsilon,kappa_tidle,sigma_alpha_tidle)


r=sparse(size(DGT,2)*N_loc,size(DGT,2)*N_loc);

%% Interior edges
for i=1:length(F_HP_i)
    n=F_HP_i(i);ele=DGE(5,n);ele_neighbor=DGE(5,DGE(6,n));
    begin_point=DGE(1:2,n);end_point=DGE(3:4,n);
    edge=sqrt((end_point(1)-begin_point(1)).^2+(end_point(2)-begin_point(2)).^2);
    
    b=zeros(N_loc,N_loc);
    b_2=zeros(N_loc,N_loc);
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)*feval(integrand_func_name_1,begin_point(1),Gauss_point_local_1D(k)...
                ,DG_edge_flag(2:4,n),DG_edge_flag(2:4,n),basis_type,eipsilon(1,ele),eipsilon(1,ele),kappa_tidle,edge,sigma_alpha_tidle);
            b_2=b_2+Gauss_coefficient_local_1D(k)*feval(integrand_func_name_1,begin_point(1),Gauss_point_local_1D(k)...
                ,DG_edge_flag(2:4,n),DG_edge_flag(2:4,DGE(6,n)),basis_type,eipsilon(1,ele),eipsilon(1,ele_neighbor),kappa_tidle,edge,sigma_alpha_tidle);
        end
    elseif begin_point(2)==end_point(2) % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)*feval(integrand_func_name_1,Gauss_point_local_1D(k),begin_point(2)...
                ,DG_edge_flag(2:4,n),DG_edge_flag(2:4,n),basis_type,eipsilon(1,ele),eipsilon(1,ele),kappa_tidle,edge,sigma_alpha_tidle);
            b_2=b_2+Gauss_coefficient_local_1D(k)*feval(integrand_func_name_1,Gauss_point_local_1D(k),begin_point(2)...
                ,DG_edge_flag(2:4,n),DG_edge_flag(2:4,DGE(6,n)),basis_type,eipsilon(1,ele),eipsilon(1,ele_neighbor),kappa_tidle,edge,sigma_alpha_tidle);
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
                ,DG_edge_flag(2:4,n),DG_edge_flag(2:4,n),basis_type,eipsilon(1,ele),eipsilon(1,ele),kappa_tidle,edge,sigma_alpha_tidle);
            b_2=b_2+Gauss_coefficient_local_1D(k)*Jacobi*feval(integrand_func_name_1,x,y...
                ,DG_edge_flag(2:4,n),DG_edge_flag(2:4,DGE(6,n)),basis_type,eipsilon(1,ele),eipsilon(1,ele_neighbor),kappa_tidle,edge,sigma_alpha_tidle);
        end
    end
    
    eipsilon_min=min(eipsilon(1,ele),eipsilon(1,ele_neighbor));
    eipsilon_max=max(eipsilon(1,ele),eipsilon(1,ele_neighbor));
    theta=eipsilon_min/eipsilon_max;

    for trial=1:N_loc
        for test=1:N_loc
            r((ele-1)*N_loc+test,(ele-1)*N_loc+trial)=r((ele-1)*N_loc+test,(ele-1)*N_loc+trial)+theta*b(trial,test);
            r((ele_neighbor-1)*N_loc+test,(ele-1)*N_loc+trial)=r((ele_neighbor-1)*N_loc+test,(ele-1)*N_loc+trial)+theta*b_2(trial,test);
        end
    end
end