function [r1,r2]=generate_load_vector_nonlinear_term_1_2(DGT,DGE,DG_edge_flag,boundary_edges,basis_type,fun_u1,fun_u2...
                                                         ,u1_h,u2_h,Gauss_coefficient_reference_1D,Gauss_point_reference_1D)

N_loc=(basis_type+1)*(basis_type+2)/2;
r1=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,1);
r2=sparse(size(DGT,2)*(basis_type+1)*(basis_type+2)/2,1);

for i=1:length(boundary_edges)
    j=boundary_edges(i);ele=DGE(5,j);
    begin_point=DGE(1:2,j);end_point=DGE(3:4,j);middle_point=(begin_point+end_point)/2;
    normal=DG_edge_flag(2:3,j);
    solution_1_ele=u1_h((ele-1)*N_loc+1:ele*N_loc,1);solution_2_ele=u2_h((ele-1)*N_loc+1:ele*N_loc,1);
    
    value=fe_solution(middle_point(1),middle_point(2),solution_1_ele,0,0,basis_type)*normal(1)...
         +fe_solution(middle_point(1),middle_point(2),solution_2_ele,0,0,basis_type)*normal(2);
     
    b_1=zeros((basis_type+1)*(basis_type+2)/2,1);
    b_2=zeros((basis_type+1)*(basis_type+2)/2,1);
    b_3=zeros((basis_type+1)*(basis_type+2)/2,1);
    b_4=zeros((basis_type+1)*(basis_type+2)/2,1);
    if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b_1=b_1+Gauss_coefficient_local_1D(k)...
                *(feval(fun_u1,begin_point(1),Gauss_point_local_1D(k))*normal(1)...
                 +feval(fun_u2,begin_point(1),Gauss_point_local_1D(k))*normal(2))...
                *feval(fun_u1,begin_point(1),Gauss_point_local_1D(k))...
                *local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,0,0)';
            b_2=b_2+Gauss_coefficient_local_1D(k)...
                *(feval(fun_u1,begin_point(1),Gauss_point_local_1D(k))*normal(1)...
                 +feval(fun_u2,begin_point(1),Gauss_point_local_1D(k))*normal(2))...
                *feval(fun_u2,begin_point(1),Gauss_point_local_1D(k))...
                *local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,0,0)';
           b_3=b_3+Gauss_coefficient_local_1D(k)...
               *(fe_solution(begin_point(1),Gauss_point_local_1D(k),solution_1_ele,0,0,basis_type)*normal(1)...
                +fe_solution(begin_point(1),Gauss_point_local_1D(k),solution_2_ele,0,0,basis_type)*normal(2))...
                *feval(fun_u1,begin_point(1),Gauss_point_local_1D(k))...
                *local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,0,0)';
           b_4=b_4+Gauss_coefficient_local_1D(k)...
               *(fe_solution(begin_point(1),Gauss_point_local_1D(k),solution_1_ele,0,0,basis_type)*normal(1)...
                +fe_solution(begin_point(1),Gauss_point_local_1D(k),solution_2_ele,0,0,basis_type)*normal(2))...
                *feval(fun_u2,begin_point(1),Gauss_point_local_1D(k))...
                *local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,0,0)'; 
        end
    else  % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b_1=b_1+Gauss_coefficient_local_1D(k)...
                *(feval(fun_u1,Gauss_point_local_1D(k),begin_point(2))*normal(1)...
                 +feval(fun_u2,Gauss_point_local_1D(k),begin_point(2))*normal(2))...
                *feval(fun_u1,Gauss_point_local_1D(k),begin_point(2))...
                *local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,0,0)';
            b_2=b_2+Gauss_coefficient_local_1D(k)...
                *(feval(fun_u1,Gauss_point_local_1D(k),begin_point(2))*normal(1)...
                 +feval(fun_u2,Gauss_point_local_1D(k),begin_point(2))*normal(2))...
                *feval(fun_u2,Gauss_point_local_1D(k),begin_point(2))...
                *local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,0,0)';
            b_3=b_3+Gauss_coefficient_local_1D(k)...
               *(fe_solution(Gauss_point_local_1D(k),begin_point(2),solution_1_ele,0,0,basis_type)*normal(1)...
                +fe_solution(Gauss_point_local_1D(k),begin_point(2),solution_2_ele,0,0,basis_type)*normal(2))...
                *feval(fun_u1,Gauss_point_local_1D(k),begin_point(2))...
                *local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,0,0)';
            b_4=b_4+Gauss_coefficient_local_1D(k)...
               *(fe_solution(Gauss_point_local_1D(k),begin_point(2),solution_1_ele,0,0,basis_type)*normal(1)...
                +fe_solution(Gauss_point_local_1D(k),begin_point(2),solution_2_ele,0,0,basis_type)*normal(2))...
                *feval(fun_u2,Gauss_point_local_1D(k),begin_point(2))...
                *local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,0,0)'; 
        end
    end
    
    for beta=1:(basis_type+1)*(basis_type+2)/2
        r1((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,1)=r1((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,1)-0.5*b_1(beta,1);
        r2((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,1)=r2((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,1)-0.5*b_2(beta,1);
        if value<0
            r1((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,1)=r1((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,1)-b_3(beta,1);
            r2((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,1)=r2((ele-1)*(basis_type+1)*(basis_type+2)/2+beta,1)-b_4(beta,1);
        end
    end
end
