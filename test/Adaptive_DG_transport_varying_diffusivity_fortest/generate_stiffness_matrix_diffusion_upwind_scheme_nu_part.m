function r=generate_stiffness_matrix_diffusion_upwind_scheme_nu_part(DGT,DGE,DG_edge_flag,basis_type,N_loc,F_HP_i,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                                    ,fun_name_2,fun_name_3,fun_name_4,eipsilon,kappa_tidle)

% Because on the interface F_HP_i, beta*n_F<0, and the normal vector is
% oriented from the element of large eipsilon to that of small eipsilon, so
% the upwind value must be taken on the element of small eipsilon, i.e.,
% the element of DG_edge_flag(4,*)=-1!
r=sparse(size(DGT,2)*N_loc,size(DGT,2)*N_loc);

%%
interior_edges=F_HP_i;

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
    j_1=edge_brother(1,i); j_2=edge_brother(2,i); 
    ele_1=DGE(5,j_1); ele_2=DGE(5,j_2);
    begin_point=DGE(1:2,j_1);end_point=DGE(3:4,j_1);
    edge=sqrt((end_point(1)-begin_point(1)).^2+(end_point(2)-begin_point(2)).^2);
    if DG_edge_flag(4,j_1)==-1
        upwind_ele=ele_1;
    else
        upwind_ele=ele_2;
    end
    
    b_1=zeros(N_loc,N_loc);b_2=zeros(N_loc,N_loc);b_3=zeros(N_loc,N_loc);b_4=zeros(N_loc,N_loc);
    b_5=zeros(N_loc,N_loc);b_6=zeros(N_loc,N_loc);b_7=zeros(N_loc,N_loc);b_8=zeros(N_loc,N_loc);
    normal_vector=DG_edge_flag(2:3,edge_brother(2,i));
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b_1=b_1+Gauss_coefficient_local_1D(k)*feval(fun_name_3,begin_point(1),Gauss_point_local_1D(k),normal_vector,basis_type,eipsilon(1,upwind_ele))*DG_edge_flag(4,j_1);
            b_2=b_2+Gauss_coefficient_local_1D(k)*feval(fun_name_3,begin_point(1),Gauss_point_local_1D(k),normal_vector,basis_type,eipsilon(1,upwind_ele))*DG_edge_flag(4,j_2);
            
            b_3=b_3+Gauss_coefficient_local_1D(k)*feval(fun_name_4,begin_point(1),Gauss_point_local_1D(k),normal_vector,basis_type,eipsilon(1,upwind_ele))*DG_edge_flag(4,j_1);
            b_4=b_4+Gauss_coefficient_local_1D(k)*feval(fun_name_4,begin_point(1),Gauss_point_local_1D(k),normal_vector,basis_type,eipsilon(1,upwind_ele))*DG_edge_flag(4,j_2);
            
            b_5=b_5+1/edge*Gauss_coefficient_local_1D(k)*feval(fun_name_2,begin_point(1),Gauss_point_local_1D(k),basis_type,0,0)*DG_edge_flag(4,j_1)*DG_edge_flag(4,j_1);
            b_6=b_6+1/edge*Gauss_coefficient_local_1D(k)*feval(fun_name_2,begin_point(1),Gauss_point_local_1D(k),basis_type,0,0)*DG_edge_flag(4,j_2)*DG_edge_flag(4,j_1);
            b_7=b_7+1/edge*Gauss_coefficient_local_1D(k)*feval(fun_name_2,begin_point(1),Gauss_point_local_1D(k),basis_type,0,0)*DG_edge_flag(4,j_1)*DG_edge_flag(4,j_2);
            b_8=b_8+1/edge*Gauss_coefficient_local_1D(k)*feval(fun_name_2,begin_point(1),Gauss_point_local_1D(k),basis_type,0,0)*DG_edge_flag(4,j_2)*DG_edge_flag(4,j_2);
        end
    elseif begin_point(2)==end_point(2) % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b_1=b_1+Gauss_coefficient_local_1D(k)*feval(fun_name_3,Gauss_point_local_1D(k),begin_point(2),normal_vector,basis_type,eipsilon(1,upwind_ele))*DG_edge_flag(4,j_1);
            b_2=b_2+Gauss_coefficient_local_1D(k)*feval(fun_name_3,Gauss_point_local_1D(k),begin_point(2),normal_vector,basis_type,eipsilon(1,upwind_ele))*DG_edge_flag(4,j_2);
            
            b_3=b_3+Gauss_coefficient_local_1D(k)*feval(fun_name_4,Gauss_point_local_1D(k),begin_point(2),normal_vector,basis_type,eipsilon(1,upwind_ele))*DG_edge_flag(4,j_1);
            b_4=b_4+Gauss_coefficient_local_1D(k)*feval(fun_name_4,Gauss_point_local_1D(k),begin_point(2),normal_vector,basis_type,eipsilon(1,upwind_ele))*DG_edge_flag(4,j_2);
            
            b_5=b_5+1/edge*Gauss_coefficient_local_1D(k)*feval(fun_name_2,Gauss_point_local_1D(k),begin_point(2),basis_type,0,0)*DG_edge_flag(4,j_1)*DG_edge_flag(4,j_1);
            b_6=b_6+1/edge*Gauss_coefficient_local_1D(k)*feval(fun_name_2,Gauss_point_local_1D(k),begin_point(2),basis_type,0,0)*DG_edge_flag(4,j_2)*DG_edge_flag(4,j_1);
            b_7=b_7+1/edge*Gauss_coefficient_local_1D(k)*feval(fun_name_2,Gauss_point_local_1D(k),begin_point(2),basis_type,0,0)*DG_edge_flag(4,j_1)*DG_edge_flag(4,j_2);
            b_8=b_8+1/edge*Gauss_coefficient_local_1D(k)*feval(fun_name_2,Gauss_point_local_1D(k),begin_point(2),basis_type,0,0)*DG_edge_flag(4,j_2)*DG_edge_flag(4,j_2);
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
            b_1=b_1+Jacobi*Gauss_coefficient_local_1D(k)*feval(fun_name_3,x,y,normal_vector,basis_type,eipsilon(1,upwind_ele))*DG_edge_flag(4,j_1);
            b_2=b_2+Jacobi*Gauss_coefficient_local_1D(k)*feval(fun_name_3,x,y,normal_vector,basis_type,eipsilon(1,upwind_ele))*DG_edge_flag(4,j_2);
            
            b_3=b_3+Jacobi*Gauss_coefficient_local_1D(k)*feval(fun_name_4,x,y,normal_vector,basis_type,eipsilon(1,upwind_ele))*DG_edge_flag(4,j_1);
            b_4=b_4+Jacobi*Gauss_coefficient_local_1D(k)*feval(fun_name_4,x,y,normal_vector,basis_type,eipsilon(1,upwind_ele))*DG_edge_flag(4,j_2);
            
            b_5=b_5+1/edge*Jacobi*Gauss_coefficient_local_1D(k)*feval(fun_name_2,x,y,basis_type,0,0)*DG_edge_flag(4,j_1)*DG_edge_flag(4,j_1);
            b_6=b_6+1/edge*Jacobi*Gauss_coefficient_local_1D(k)*feval(fun_name_2,x,y,basis_type,0,0)*DG_edge_flag(4,j_2)*DG_edge_flag(4,j_1);
            b_7=b_7+1/edge*Jacobi*Gauss_coefficient_local_1D(k)*feval(fun_name_2,x,y,basis_type,0,0)*DG_edge_flag(4,j_1)*DG_edge_flag(4,j_2);
            b_8=b_8+1/edge*Jacobi*Gauss_coefficient_local_1D(k)*feval(fun_name_2,x,y,basis_type,0,0)*DG_edge_flag(4,j_2)*DG_edge_flag(4,j_2);
        end
    end
    
    eipsilon_min=min(eipsilon(1,ele_1),eipsilon(1,ele_2));
    eipsilon_max=max(eipsilon(1,ele_1),eipsilon(1,ele_2));
    theta=eipsilon_min/eipsilon_max;
    
    sigma_nu_tidle=eipsilon_min;
    for trial=1:N_loc
        for test=1:N_loc
            r((ele_1-1)*N_loc+test,(upwind_ele-1)*N_loc+trial)=r((ele_1-1)*N_loc+test,(upwind_ele-1)*N_loc+trial)+(1-theta)*b_1(trial,test);
            r((ele_2-1)*N_loc+test,(upwind_ele-1)*N_loc+trial)=r((ele_2-1)*N_loc+test,(upwind_ele-1)*N_loc+trial)+(1-theta)*b_2(trial,test);
            
            r((upwind_ele-1)*N_loc+test,(ele_1-1)*N_loc+trial)=r((upwind_ele-1)*N_loc+test,(ele_1-1)*N_loc+trial)+kappa_tidle*(1-theta)*b_3(trial,test);
            r((upwind_ele-1)*N_loc+test,(ele_2-1)*N_loc+trial)=r((upwind_ele-1)*N_loc+test,(ele_2-1)*N_loc+trial)+kappa_tidle*(1-theta)*b_4(trial,test);
            
            r((ele_1-1)*N_loc+test,(ele_1-1)*N_loc+trial)=r((ele_1-1)*N_loc+test,(ele_1-1)*N_loc+trial)+sigma_nu_tidle*(1-theta)*b_5(trial,test);
            r((ele_2-1)*N_loc+test,(ele_1-1)*N_loc+trial)=r((ele_2-1)*N_loc+test,(ele_1-1)*N_loc+trial)+sigma_nu_tidle*(1-theta)*b_6(trial,test);
            r((ele_1-1)*N_loc+test,(ele_2-1)*N_loc+trial)=r((ele_1-1)*N_loc+test,(ele_2-1)*N_loc+trial)+sigma_nu_tidle*(1-theta)*b_7(trial,test);
            r((ele_2-1)*N_loc+test,(ele_2-1)*N_loc+trial)=r((ele_2-1)*N_loc+test,(ele_2-1)*N_loc+trial)+sigma_nu_tidle*(1-theta)*b_8(trial,test);
        end
    end
end


    
    
    