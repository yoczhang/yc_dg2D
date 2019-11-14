function r=integrand_3_on_edges(x,y,normal_vector,basis_type,eipsilon)
%The integrand of the integral on the inner edges


r=-eipsilon*(normal_vector(1)*triangular_local_basis(x,y,basis_type,1,0)...
                      +normal_vector(2)*triangular_local_basis(x,y,basis_type,0,1))'...
                     *triangular_local_basis(x,y,basis_type,0,0);