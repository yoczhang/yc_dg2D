function r=generate_stiffness_matrix_from_convection_term_1(DGE,DGT,DG_edge_flag,basis_type,interior_edges,flag,vector_u_1,vector_u_2,integrand_func_name...
                                                           ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D)

% Because we know the orientation of the velocity (-0.2,-0.1), the upwind value of
% concentration must be taken as the limitation of that on the larger element

r=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,size(DGT,2)*(basis_type+1)*(basis_type+2)/2);

%%
% The 1st line of edge_brother is the small element
% its 2nd line is the large element
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

%% 
for i=1:size(edge_brother,2)
    ele_large=DGE(5,edge_brother(2,i));ele_small=DGE(5,edge_brother(1,i));

    begin_point=DGE(1:2,edge_brother(1,i));end_point=DGE(3:4,edge_brother(1,i));
    
    b=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    b_2=zeros((basis_type+1)*(basis_type+2)/2,(basis_type+1)*(basis_type+2)/2);
    if flag(ele_large)/flag(ele_small)==1
        if flag(ele_large)==1
            vector_u=vector_u_1;
        else
            vector_u=vector_u_2;
        end
    else
        vector_u=(vector_u_1+vector_u_2)/2;
    end
    
    temp=vector_u(1)*DG_edge_flag(2,edge_brother(2,i))+vector_u(2)*DG_edge_flag(3,edge_brother(2,i));
    if temp<0
        ele_trial=ele_small;temp_1=-1;
        ele_test=ele_large;temp_2=1;
    else
        ele_trial=ele_large;temp_1=1;
        ele_test=ele_small;temp_2=-1;
    end
    
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+temp_1*temp*Gauss_coefficient_local_1D(k)*feval(integrand_func_name,begin_point(1),Gauss_point_local_1D(k),basis_type);
            b_2=b_2+temp_2*temp*Gauss_coefficient_local_1D(k)*feval(integrand_func_name,begin_point(1),Gauss_point_local_1D(k),basis_type);
        end
    elseif begin_point(2)==end_point(2) % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+temp_1*temp*Gauss_coefficient_local_1D(k)*feval(integrand_func_name,Gauss_point_local_1D(k),begin_point(2),basis_type);
            b_2=b_2+temp_2*temp*Gauss_coefficient_local_1D(k)*feval(integrand_func_name,Gauss_point_local_1D(k),begin_point(2),basis_type);
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
            b=b+Jacobi*temp_1*temp*Gauss_coefficient_local_1D(k)*feval(integrand_func_name,x,y,basis_type);
            b_2=b_2+Jacobi*temp_2*temp*Gauss_coefficient_local_1D(k)*feval(integrand_func_name,x,y,basis_type);
        end
    end
    
    for alpha=1:(basis_type+1)*(basis_type+2)/2
        for beta=1:(basis_type+1)*(basis_type+2)/2
            r((ele_trial-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele_trial-1)*(basis_type+1)*(basis_type+2)/2+alpha)...
                =r((ele_trial-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele_trial-1)*(basis_type+1)*(basis_type+2)/2+alpha)+b(alpha,beta);
            r((ele_test-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele_trial-1)*(basis_type+1)*(basis_type+2)/2+alpha)...
                =r((ele_test-1)*(basis_type+1)*(basis_type+2)/2+beta,(ele_trial-1)*(basis_type+1)*(basis_type+2)/2+alpha)+b_2(alpha,beta);
        end
    end
end

