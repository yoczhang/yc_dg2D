function r=generate_load_vector_inflow_boundary_edges(DGT,DGE,basis_type,N_loc,F_in,Gauss_coefficient_reference_1D,Gauss_point_reference_1D)

r=sparse(size(DGT,2)*N_loc,1);

for i=1:length(F_in)
    j=F_in(i);
    begin_point=DGE(1:2,j);end_point=DGE(3:4,j);
    ele=DGE(5,j);
    
    b=zeros(N_loc,1);
   if begin_point(1)==end_point(1) %vertical edge
        lower_bound=min(begin_point(2),end_point(2));
        upper_bound=max(begin_point(2),end_point(2));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)*triangular_local_basis(begin_point(1),Gauss_point_local_1D(k),basis_type,0,0)';
        end
   else  % begin_point(2)==end_point(2) % horizontal edge
        lower_bound=min(begin_point(1),end_point(1));
        upper_bound=max(begin_point(1),end_point(1));
        [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
        for k=1:length(Gauss_coefficient_local_1D)
            b=b+Gauss_coefficient_local_1D(k)*triangular_local_basis(Gauss_point_local_1D(k),begin_point(2),basis_type,0,0)';
        end
   end
   
   for test=1:N_loc
       r((ele-1)*N_loc+test,1)=r((ele-1)*N_loc+test,1)+b(test,1);
    end
end

