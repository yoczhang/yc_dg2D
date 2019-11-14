function r=generate_stiffness_matrix_nonlinear_term_3_4(DGT,DGE,DG_edge_flag,fun_name_2,fun_name_3,fun_name_4,fun_name_5,basis_type,u1_h,u2_h...
                                                       ,interior_edges,boundary_edges,Gauss_coefficient_reference_1D,Gauss_point_reference_1D)

N_loc=(basis_type+1)*(basis_type+2)/2;
r=sparse(size(DGT,2)*N_loc,size(DGT,2)*N_loc);

%%
% The 1st line of edge_brother is the small element
% its 2nd line is the large element
edge_brother=zeros(2,length(interior_edges)/2);
j=1;
for i=1:length(interior_edges)
    if ismember(interior_edges(i),edge_brother)==0
        edge_brother(1,j)=interior_edges(i);edge_brother(2,j)=DGE(6,interior_edges(i));
    else
        continue
    end  
    j=j+1;
end
clear j

%% Interior edges
for i=1:size(edge_brother,2)
    ele_large=DGE(5,edge_brother(2,i));ele_small=DGE(5,edge_brother(1,i));
    begin_point=DGE(1:2,edge_brother(1,i));end_point=DGE(3:4,edge_brother(1,i));middle_point=(begin_point+end_point)/2;
    normal=DG_edge_flag(2:3,edge_brother(1,i));
    solution_1_large_ele=u1_h((ele_large-1)*N_loc+1:ele_large*N_loc,1);solution_2_large_ele=u2_h((ele_large-1)*N_loc+1:ele_large*N_loc,1);
    solution_1_small_ele=u1_h((ele_small-1)*N_loc+1:ele_small*N_loc,1);solution_2_small_ele=u2_h((ele_small-1)*N_loc+1:ele_small*N_loc,1);
    
    value=0.5*((fe_solution(middle_point(1),middle_point(2),solution_1_small_ele,0,0,basis_type)...
               +fe_solution(middle_point(1),middle_point(2),solution_1_large_ele,0,0,basis_type))*normal(1)...
              +(fe_solution(middle_point(1),middle_point(2),solution_2_small_ele,0,0,basis_type)...
               +fe_solution(middle_point(1),middle_point(2),solution_2_large_ele,0,0,basis_type))*normal(2));
    if value<0
        test_ele=ele_large;
    else
        test_ele=ele_small;
    end
    
    b_1=zeros(N_loc,N_loc);
    b_2=zeros(N_loc,N_loc);
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b_1=b_1+Gauss_coefficient_local_1D(k)...
                   *feval(fun_name_2,begin_point(1),Gauss_point_local_1D(k),solution_1_large_ele,solution_1_small_ele,solution_2_large_ele,solution_2_small_ele,basis_type,normal);
            b_2=b_2+Gauss_coefficient_local_1D(k)...
                   *feval(fun_name_3,begin_point(1),Gauss_point_local_1D(k),solution_1_large_ele,solution_1_small_ele,solution_2_large_ele,solution_2_small_ele,basis_type,normal);
        end
    elseif begin_point(2)==end_point(2) % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b_1=b_1+Gauss_coefficient_local_1D(k)...
                   *feval(fun_name_2,Gauss_point_local_1D(k),begin_point(2),solution_1_large_ele,solution_1_small_ele,solution_2_large_ele,solution_2_small_ele,basis_type,normal);
            b_2=b_2+Gauss_coefficient_local_1D(k)...
                   *feval(fun_name_3,Gauss_point_local_1D(k),begin_point(2),solution_1_large_ele,solution_1_small_ele,solution_2_large_ele,solution_2_small_ele,basis_type,normal);
        end
    else
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        slope=(end_point(2)-begin_point(2))/(end_point(1)-begin_point(1));
        Jacobi=sqrt(1+slope^2);
        for k=1:length(Gauss_coefficient_local_1D)
            x=Gauss_point_local_1D(k);y=slope*(x-begin_point(1))+begin_point(2);
            b_1=b_1+Gauss_coefficient_local_1D(k)*Jacobi...
                *feval(fun_name_2,x,y,solution_1_large_ele,solution_1_small_ele,solution_2_large_ele,solution_2_small_ele,basis_type,normal);
            b_2=b_2+Gauss_coefficient_local_1D(k)*Jacobi...
                *feval(fun_name_3,x,y,solution_1_large_ele,solution_1_small_ele,solution_2_large_ele,solution_2_small_ele,basis_type,normal);
        end
    end
    
    for alpha=1:N_loc
        for beta=1:N_loc
            r((ele_large-1)*N_loc+beta,(ele_large-1)*N_loc+alpha)=r((ele_large-1)*N_loc+beta,(ele_large-1)*N_loc+alpha)+b_1(alpha,beta);
            r((ele_small-1)*N_loc+beta,(ele_small-1)*N_loc+alpha)=r((ele_small-1)*N_loc+beta,(ele_small-1)*N_loc+alpha)+b_1(alpha,beta);
            r((test_ele-1)*N_loc+beta,(ele_small-1)*N_loc+alpha)=r((test_ele-1)*N_loc+beta,(ele_small-1)*N_loc+alpha)+b_2(alpha,beta);
            r((test_ele-1)*N_loc+beta,(ele_large-1)*N_loc+alpha)=r((test_ele-1)*N_loc+beta,(ele_large-1)*N_loc+alpha)-b_2(alpha,beta);
        end
    end
end
clear i ele_large ele_small begin_point end_point normal value test_ele b_1 b_2

%% Boundary edges
for i=1:length(boundary_edges)
    j=boundary_edges(i);ele=DGE(5,j);
    begin_point=DGE(1:2,j);end_point=DGE(3:4,j);middle_point=(begin_point+end_point)/2;
    normal=DG_edge_flag(2:3,j);
    solution_1_ele=u1_h((ele-1)*N_loc+1:ele*N_loc,1);solution_2_ele=u2_h((ele-1)*N_loc+1:ele*N_loc,1);
    
    value=fe_solution(middle_point(1),middle_point(2),solution_1_ele,0,0,basis_type)*normal(1)...
         +fe_solution(middle_point(1),middle_point(2),solution_2_ele,0,0,basis_type)*normal(2);           
    b_1=zeros(N_loc,N_loc);
    b_2=zeros(N_loc,N_loc);
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b_1=b_1+Gauss_coefficient_local_1D(k)...
                   *feval(fun_name_4,begin_point(1),Gauss_point_local_1D(k),solution_1_ele,solution_2_ele,basis_type,normal);
            b_2=b_2+Gauss_coefficient_local_1D(k)...
                   *feval(fun_name_5,begin_point(1),Gauss_point_local_1D(k),solution_1_ele,solution_2_ele,basis_type,normal);
        end
    else % begin_point(2)==end_point(2) % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b_1=b_1+Gauss_coefficient_local_1D(k)...
                   *feval(fun_name_4,Gauss_point_local_1D(k),begin_point(2),solution_1_ele,solution_2_ele,basis_type,normal);
            b_2=b_2+Gauss_coefficient_local_1D(k)...
                   *feval(fun_name_5,Gauss_point_local_1D(k),begin_point(2),solution_1_ele,solution_2_ele,basis_type,normal);
        end
    end
    
    for alpha=1:N_loc
        for beta=1:N_loc
            r((ele-1)*N_loc+beta,(ele-1)*N_loc+alpha)=r((ele-1)*N_loc+beta,(ele-1)*N_loc+alpha)+b_1(alpha,beta);
            if middle_point(1)==0
                r((ele-1)*N_loc+beta,(ele-1)*N_loc+alpha)=r((ele-1)*N_loc+beta,(ele-1)*N_loc+alpha)+b_2(alpha,beta);
            end
        end
    end
end    
    
    
    
    