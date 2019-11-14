function r=generate_stiffness_matrix_from_convection_term(DGE,DGT,DG_edge_flag,basis_type,N_loc,F_i,F_HP_i,F_out,integrand_func_name...
                                                           ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D)
% The wind direction beta=[1, 0]

r=sparse(size(DGT,2)*N_loc,size(DGT,2)*N_loc);

%%
% The 1st line of edge_brother is the small element
% its 2nd line is the large element
interior_edges=[F_i F_HP_i];

edge_brother=zeros(2,length(interior_edges)/2);
j=1;
for i=1:length(interior_edges)
    if ismember(interior_edges(i),edge_brother)==0
        edge_brother(1,j)=interior_edges(i);
        edge_brother(2,j)=DGE(6,interior_edges(i));
    else
        continue
    end  
    j=j+1;
end
clear j

%%  Interior edges
for i=1:size(edge_brother,2)
    begin_point=DGE(1:2,edge_brother(1,i));end_point=DGE(3:4,edge_brother(1,i));
    
    b=zeros(N_loc,N_loc);b_2=zeros(N_loc,N_loc);
    
    if DG_edge_flag(4,edge_brother(2,i))==1
        active=edge_brother(2,i);passive=edge_brother(1,i);
    else
        active=edge_brother(1,i);passive=edge_brother(2,i);
    end

    if DG_edge_flag(2,edge_brother(2,i))<0  % means that beta*n_F<0
        ele_trial=DGE(5,passive);
    else
        ele_trial=DGE(5,active);
    end
    
    temp=DG_edge_flag(2,edge_brother(2,i));
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+temp*DG_edge_flag(4,edge_brother(1,i))*Gauss_coefficient_local_1D(k)*feval(integrand_func_name,begin_point(1),Gauss_point_local_1D(k),basis_type,0,0);
            b_2=b_2+temp*DG_edge_flag(4,edge_brother(2,i))*Gauss_coefficient_local_1D(k)*feval(integrand_func_name,begin_point(1),Gauss_point_local_1D(k),basis_type,0,0);
        end
    elseif begin_point(2)==end_point(2) % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+temp*DG_edge_flag(4,edge_brother(1,i))*Gauss_coefficient_local_1D(k)*feval(integrand_func_name,Gauss_point_local_1D(k),begin_point(2),basis_type,0,0);
            b_2=b_2+temp*DG_edge_flag(4,edge_brother(2,i))*Gauss_coefficient_local_1D(k)*feval(integrand_func_name,Gauss_point_local_1D(k),begin_point(2),basis_type,0,0);
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
            b=b+temp*DG_edge_flag(4,edge_brother(1,i))*Gauss_coefficient_local_1D(k)*Jacobi*feval(integrand_func_name,x,y,basis_type,0,0);
            b_2=b_2+temp*DG_edge_flag(4,edge_brother(2,i))*Gauss_coefficient_local_1D(k)*Jacobi*feval(integrand_func_name,x,y,basis_type,0,0);
        end
    end
    
    for trial=1:N_loc
        for test=1:N_loc
            r((DGE(5,edge_brother(1,i))-1)*N_loc+test,(ele_trial-1)*N_loc+trial)...
                =r((DGE(5,edge_brother(1,i))-1)*N_loc+test,(ele_trial-1)*N_loc+trial)+b(trial,test);
            r((DGE(5,edge_brother(2,i))-1)*N_loc+test,(ele_trial-1)*N_loc+trial)...
                =r((DGE(5,edge_brother(2,i))-1)*N_loc+test,(ele_trial-1)*N_loc+trial)+b_2(trial,test);
        end
    end
end

%% F_out
for i=1:length(F_out)
    j=F_out(i); ele=DGE(5,j);
    begin_point=DGE(1:2,j);end_point=DGE(3:4,j);
    
    b=zeros(N_loc,N_loc);
    
    temp=DG_edge_flag(2,j);
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+temp*Gauss_coefficient_local_1D(k)*feval(integrand_func_name,begin_point(1),Gauss_point_local_1D(k),basis_type,0,0);
        end
    else % begin_point(2)==end_point(2)  horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+temp*Gauss_coefficient_local_1D(k)*feval(integrand_func_name,Gauss_point_local_1D(k),begin_point(2),basis_type,0,0);
        end
    end
    
    for trial=1:N_loc
        for test=1:N_loc
            r((ele-1)*N_loc+test,(ele-1)*N_loc+trial)=r((ele-1)*N_loc+test,(ele-1)*N_loc+trial)+b(trial,test);
        end
    end
end
