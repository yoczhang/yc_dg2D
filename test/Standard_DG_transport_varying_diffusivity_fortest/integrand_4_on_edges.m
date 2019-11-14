function r=integrand_4_on_edges(x,y,normal_vector,basis_type,eipsilon)
%The integrand of the integral on the inner edges


r=triangular_local_basis(x,y,basis_type,0,0)'...
    *eipsilon*(normal_vector(1)*triangular_local_basis(x,y,basis_type,1,0)...
                      +normal_vector(2)*triangular_local_basis(x,y,basis_type,0,1));